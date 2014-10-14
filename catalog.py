import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import fitsio
import re
import os

## copying /home/howdiedoo/boss/mksamplecatslatestdr11/fitstotxtfkpNEWSYS.py 
def catfitstotxt(zmin,zmax,fdir='./',sampletag='cmass',runtag='dr12v4',cattypetag = 'Reid',NorS = 2,
                 cpopt=2,sysopt=1,fkpopt=0,DorR=0,ftot=None):
  """
  Convert data or random fits catalog (output from mksample) 
  to txt file to be read into my pair counting code.
  Given sampletag, runtag, cattypetag, NorS (2=both), will guess filenames for data and random.
  Input a filename as ftot if you have a different pattern/just want 
  to do one file.
  fdir is the directory containing the fits files.
  GOING TO MAKE THIS FXN NOT BACKWARDS COMPATIBLE!!
  I don't *think* anything regular calls this function, so it *should be ok.
  sysopt = 0 (none) or 1 (use 'WEIGHT_SYSTOT')
  fkpopt = 0 (none) or 1 (use 'WEIGHT_FKP')
  cpopt = 0 (all items in catalog equal; appropriate for 'target' catalogs which have close pair weights but we don't want to use them)
  cpopt = 1 (dd['WEIGHT_NOZ'])  [appropriate for 'ang' catalogs]
  cpopt = 2 (dd['WEIGHT_NOZ'][i] + dd['WEIGHT_CP'][i] - 1.0) [appropriate for NN upweighting scheme catalogs]
  """

  assert sysopt == 0 or sysopt == 1
  assert fkpopt == 0 or fkpopt == 1
  assert cpopt >= 0 and cpopt <=2
  assert DorR == 0 or DorR == 1

  if sampletag == 'cmass':
    assert sysopt == 1
  else:
    assert sysopt == 0

  if DorR == 1:
    assert cpopt == 0

  ## enforce sanity on catalog type and weighting scheme to protect myself from screw-ups!
  targcat = 0
  if re.search('targ',cattypetag):
    targcat = 1  # no redshifts!
    assert fkpopt == 0 
    assert cpopt == 0

  elif re.search('ang',cattypetag):
    assert fkpopt == 0 
    if DorR == 0:
      assert cpopt == 1

  else:
    if DorR == 0:
      assert cpopt == 2

    

  catappend = ''
  if sysopt == 1:
    catappend = catappend + '-wsys'
  if fkpopt == 1:
    catappend = catappend + '-wfkp'

  if ftot is not None:
    flist = [ftot]
    DorRlist = [DorR]
    syslist = [sysopt]
  ## do full pattern of data and random, N and S as indicated.
  else:
    flist = []
    foutlist = []
    DorRlist = []
    if NorS == 0:
      NSlist = ['N']
    elif NorS == 1:
      NSlist = ['S']
    else:
      NSlist = ['N','S']

    for NS in NSlist:
      for DRval in [0,1]:
        if DRval == 0:
          fname = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + '.dat.fits'
          fout = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + '.dat.txt'
        if DRval == 1:
          fname = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + '.ran.fits'
          fout = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + '.ran.txt'
        flist.append(fdir + fname)
        foutlist.append(fdir + fout)
        DorRlist.append(DRval)


  for ff, fout, DRval in zip(flist,foutlist,DorRlist):
    dd = fitsio.read(ff)
    print 'hi',ff,fout,DRval
    ## logic for weighting data/randoms.
    if cpopt == 0 or DRval == 1:
      wgtvec = np.zeros(len(dd)) + 1.
    if cpopt == 1 and DRval == 0:
      wgtvec = dd['WEIGHT_NOZ']
    if cpopt == 2 and DRval == 0:
      wgtvec = (dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0) 
    if sysopt == 1 and DRval == 0:
      wgtvec = wgtvec * dd['WEIGHT_SYSTOT']
    if fkpopt == 1:
      wgtvec = wgtvec * dd['WEIGHT_FKP']
    if DRval == 0:
      h, ii = np.histogram(dd['IMATCH'][:],bins=np.arange(-0.5,10.5,1))
      if targcat == 0:
        assert h[0] == 0
        assert (h[3:] == 0).all()
      assert (np.fabs(dd['WEIGHT_STAR']*dd['WEIGHT_SEEING'] - dd['WEIGHT_SYSTOT'])/dd['WEIGHT_SYSTOT'] < 2.0e-5).all()
      if sysopt == 1:
        syswgt = dd['WEIGHT_SYSTOT']*(dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0)
    ofp = open(fout,'w')
    if targcat == 0:
      for i in range(len(dd)):
        if(dd['Z'][i] < zmin or dd['Z'][i] > zmax):
          continue
        ofp.write('%.12e %.12e %.6e %.6e\n' % (dd['RA'][i],dd['DEC'][i],dd['Z'][i],wgtvec[i])) 
    else: ## do target catalog.
      for i in range(len(dd)):
        ofp.write('%.12e %.12e %.6e\n' % (dd['RA'][i],dd['DEC'][i],wgtvec[i])) 

    ofp.close()
    print 'finished printing ',ff,'with this many elements:',len(dd), 'sysopt = ',sysopt, 'fkpopt = ',fkpopt


def writexiparamfile(pfname,runp,DRopt,Dftype=1,Rftype=1):
  """
  runp is a dictionary with all hte relevant parameters.
  targcat = 1 for catalog of targets, 0 otherwise.  targets have different file format (3), usual ra,dec,z,wgt (1).
  """

  mykeylist = ['omfid','hfid','binfname','zmin','zmax','foutbase','Dfilename','Rfilename']
  #defaultlist = [0.292, 0.69, "xibinfiles/bin_xismu1Mpcohtest.txt", 0.43, 0.7, "outputdr12/testing","tmp",'tmp']

  ofp = open(pfname,'w')
  ## things we're holding fixed for now.
  ofp.write('DRopt = %d\n' % DRopt)
  ofp.write('radeczorsim = 0\n')
  ofp.write('unitsMpc = 0\n')
  ofp.write('Dftype = %d\n' % (Dftype))
  ofp.write('Rftype = %d\n' % (Rftype))
  #how are zmin/zmax assigned in xi?  need to set them to some big range in the params file.

  for kk in mykeylist:
    try:
      aa = runp[kk]
    except:
      #if kk is not 'omfid' and kk is not 'hfid':
      print 'key not found, aborting!'
      sys.exit(1)

    if kk == 'ndownRR' or kk == 'rngenseed' or kk == 'ndownDD': continue # deal with these later.

    if type(runp[kk]) is str:
      ofp.write('%s = %s\n' % (kk,runp[kk]))
    if type(runp[kk]) is int:
      ofp.write('%s = %d\n' % (kk,runp[kk]))
    if type(runp[kk]) is float:
      ofp.write('%s = %e\n' % (kk,runp[kk]))

  myflag = 0

  if (DRopt == 3 or DRopt == 12 or DRopt == 14)  and runp['ndownRR'] - 1.0 > 0.001:
    ofp.write('%s = %f\n' % ('ndownRR',runp['ndownRR']))
    myflag = 1
  if (DRopt == 13 or DRopt == 14)  and runp['ndownDD'] - 1.0 > 0.001:
    ofp.write('%s = %f\n' % ('ndownDD',runp['ndownDD']))
    myflag = 1
  if myflag == 1:  # write a random seed if we're going to downsample.
    if runp['rngenseed'] > 0:
      ofp.write('%s = %d\n' % ('rngenseed',runp['rngenseed']))
    else:
      myrand = np.random.randint(1,2**31-1,1)
      assert myrand > 0
      ofp.write('%s = %d\n' % ('rngenseed',myrand))
  ofp.close()

def catfitstotxtpatch(zmin,zmax,fdir='./',sampletag='cmass',runtag='dr12v4',cattypetag = 'Reid',NorS = 2,
                 cpopt=2,sysopt=1,fkpopt=0,DorR=0,ftot=None,ramin=-720.,ramax=720.,decmin=-1000.,decmax=1000.):
  """
  copied catfitstotxt but added cut to a specified ra,dec box for testing.
  """

  assert sysopt == 0 or sysopt == 1
  assert fkpopt == 0 or fkpopt == 1
  assert cpopt >= 0 and cpopt <=2
  assert DorR == 0 or DorR == 1

  if sampletag == 'cmass':
    assert sysopt == 1
  else:
    assert sysopt == 0

  if DorR == 1:
    assert cpopt == 0

  ## enforce sanity on catalog type and weighting scheme to protect myself from screw-ups!
  targcat = 0
  if re.search('targ',cattypetag):
    targcat = 1  # no redshifts!
    assert fkpopt == 0 
    assert cpopt == 0

  elif re.search('ang',cattypetag):
    assert fkpopt == 0 
    if DorR == 0:
      assert cpopt == 1

  else:
    if DorR == 0:
      assert cpopt == 2

    

  catappend = ''
  if sysopt == 1:
    catappend = catappend + '-wsys'
  if fkpopt == 1:
    catappend = catappend + '-wfkp'

  catappend = catappend + '-patchtest'

  if ftot is not None:
    flist = [ftot]
    DorRlist = [DorR]
    syslist = [sysopt]
  ## do full pattern of data and random, N and S as indicated.
  else:
    flist = []
    foutlist = []
    DorRlist = []
    if NorS == 0:
      NSlist = ['N']
    elif NorS == 1:
      NSlist = ['S']
    else:
      NSlist = ['N','S']

    for NS in NSlist:
      for DRval in [0,1]:
        if DRval == 0:
          fname = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + '.dat.fits'
          fout = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + '.dat.txt'
        if DRval == 1:
          fname = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + '.ran.fits'
          fout = sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + '.ran.txt'
        flist.append(fdir + fname)
        foutlist.append(fdir + fout)
        DorRlist.append(DRval)


  for ff, fout, DRval in zip(flist,foutlist,DorRlist):
    dd = fitsio.read(ff)
    print 'hi',ff,fout,DRval
    ## logic for weighting data/randoms.
    if cpopt == 0 or DRval == 1:
      wgtvec = np.zeros(len(dd)) + 1.
    if cpopt == 1 and DRval == 0:
      wgtvec = dd['WEIGHT_NOZ']
    if cpopt == 2 and DRval == 0:
      wgtvec = (dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0) 
    if sysopt == 1 and DRval == 0:
      wgtvec = wgtvec * dd['WEIGHT_SYSTOT']
    if fkpopt == 1:
      wgtvec = wgtvec * dd['WEIGHT_FKP']
    if DRval == 0:
      h, ii = np.histogram(dd['IMATCH'][:],bins=np.arange(-0.5,10.5,1))
      if targcat == 0:
        assert h[0] == 0
        assert (h[3:] == 0).all()
      assert (np.fabs(dd['WEIGHT_STAR']*dd['WEIGHT_SEEING'] - dd['WEIGHT_SYSTOT'])/dd['WEIGHT_SYSTOT'] < 2.0e-5).all()
      if sysopt == 1:
        syswgt = dd['WEIGHT_SYSTOT']*(dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0)
    ofp = open(fout,'w')
    if targcat == 0:
      for i in range(len(dd)):
        if(dd['Z'][i] < zmin or dd['Z'][i] > zmax):
          continue
        if(dd['RA'][i] < ramin or dd['RA'][i] > ramax or dd['DEC'][i] < decmin or dd['DEC'][i] > decmax):
          continue

        ofp.write('%.12e %.12e %.6e %.6e\n' % (dd['RA'][i],dd['DEC'][i],dd['Z'][i],wgtvec[i])) 
    else: ## do target catalog.
      for i in range(len(dd)):
        if(dd['RA'][i] < ramin or dd['RA'][i] > ramax or dd['DEC'][i] < decmin or dd['DEC'][i] > decmax):
          continue
        ofp.write('%.12e %.12e %.6e\n' % (dd['RA'][i],dd['DEC'][i],wgtvec[i])) 

    ofp.close()
    print 'finished printing ',ff,'with this many elements:',len(dd), 'sysopt = ',sysopt, 'fkpopt = ',fkpopt

def callxi(zmin,zmax,ndownRR=1.0,ndownDD=1.0,rngenseed=-1,ddir="/home/howdiedoo/boss/mksamplecatslatestdr12/",sampletag='cmass',runtag='dr12v4',cattypetag='Reid',NorS = 2, sysopt=1,fkpopt=0,omfid=0.292, hfid = 0.69, whichtask = 0, runopt = 0, outdirbase = 'outputdr12',fname1pw = None, binfname = None, ztag = '',cattypetag2=None,Nsub=-1,Nsubdir = None):
  """
  whichtask = 0: xiell
  whichtask = 1: wp (compute xi(rp,rpi))
  whichtask = 2: wtheta
  whichtask = 3: Hogg spec-im cross-correlation.
  runopt = 0: just print param file for calling elsewhere.
  runopt = 1: actually call xi and produce hte result.
  runopt = 2: for DRopt = 3 only, submit with nohup to run in background.
  runopt = 3: for DRopt = 1 and 2 only, run the command.
  Assuming that you are running this code from a directory containing the following folders:
  outputdirbase + '-xiell', '-xigrid', '-wtheta', '-wpcross'
  xibinfiles with all the appropriate bin files.
  Also assuming the ./xi executable is in this folder.
  ndownRR and rngenseed only written if DRopt = 3 or DRopt == 12,14 and rngenseed > 0 and (ndownRR - 1.0) > 0.001
  ndownDD only written if DRopt = 13,14
  cattypetag2 specifies tag for imaging catalog for the Hogg method to get wp.
  """

  bootopt = 0
  if Nsub > 0:
    bootopt = 1
    assert Nsubdir is not None
  print 'hi beth bootopt',bootopt

  assert zmin < zmax ## this should work for ang case, doesn't matter what they are.

  targcat = 0
  Dftype = 1
  Rftype = 1
  if re.search('targ',cattypetag):
    targcat = 1  # no redshifts!
    Dftype = 3
  if cattypetag2 is not None:
    if re.search('targ',cattypetag2):
      Rftype = 3

  cmdlist = []

  catappend = ''
  if sysopt == 1:
    catappend = catappend + '-wsys'
  if fkpopt == 1:
    catappend = catappend + '-wfkp'

  catappend = catappend

  runp = {}
  runp['omfid'] = omfid
  runp['hfid'] = hfid
  runp['zmin'] = zmin
  runp['zmax'] = zmax
  runp['ndownRR'] = ndownRR
  runp['ndownDD'] = ndownDD
  runp['rngenseed'] = rngenseed

  if whichtask == 0:
    outdir = outdirbase + '-xiell'
    runp['binfname'] = "xibinfiles/bin1xiellsmallscale.txt"
  if whichtask == 1:
    outdir = outdirbase + '-xigrid'
    runp['binfname'] = "xibinfiles/bin1wp.txt"
  if whichtask == 2:
    outdir = outdirbase + '-wtheta'
    runp['binfname'] = "xibinfiles/bin1ang.txt"
  if whichtask == 3:
    outdir = outdirbase + '-wpcross'
    runp['binfname'] = "xibinfiles/bin1wpsmall.txt"

  if NorS == 0:
    NSlist = ['N']
  elif NorS == 1:
    NSlist = ['S']
  else:
    NSlist = ['N','S']

  ## overload NSlist with all the subregion appendices
  if bootopt == 1:
    NSlist = np.arange(0,Nsub,1,dtype='int')
    mycmd = 'mkdir %s' % (outdir+'/'+Nsubdir+'/')
    print 'about to run',mycmd
    os.system(mycmd)

  DRoptlist = [1,2,3]
  if whichtask == 3:
    DRoptlist = [11,12,13,14]

  for NS in NSlist:

    ## these are for whichtask <= 2, not Hogg.
    if bootopt == 0:
      runp['Dfilename'] =  ddir + sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + '.dat.txt'
      runp['Rfilename'] =  ddir + sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + '.ran.txt'
    else:
      runp['Dfilename'] =  ddir + Nsubdir + '/' + sampletag + '-' + runtag + '-' + cattypetag + catappend + '.dat.txt' + '.%04d' % NS
      runp['Rfilename'] =  ddir + Nsubdir + '/' + sampletag + '-' + runtag + '-' + cattypetag + catappend + '.ran.txt' + '.%04d' % NS
       

    if binfname is not None:  #override defaults!
      runp['binfname'] = binfname

    ## is this generic enough NOT to overwrite myself?  I think so, cattypetag should do most of the discrimination.
    if bootopt == 0:
      runp['foutbase'] = outdir + '/' + sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend + ztag
    else:
      runp['foutbase'] = outdir + '/' +  Nsubdir + '/' + sampletag + '-' + runtag + '-' + cattypetag + catappend + ztag + '.%04d' % NS
    for DRopt in DRoptlist:
      if whichtask == 3:
        if bootopt == 0:
          dbase =  ddir + sampletag + '-' + runtag + '-' + NS + '-' + cattypetag + catappend
          rbase = ddir + sampletag + '-' + runtag + '-' + NS + '-' + cattypetag2 + catappend
          fend = ''
        else:
          dbase =  ddir + Nsubdir + '/' + sampletag + '-' + runtag +  '-' + cattypetag + catappend
          rbase = ddir + Nsubdir + '/' + sampletag + '-' + runtag  + '-' + cattypetag2 + catappend
          fend = '.%04d' % NS

        if DRopt == 11:
          runp['Dfilename'] = dbase + '.dat.txt' + fend 
          runp['Rfilename'] = rbase + '.dat.txt' + fend
        if DRopt == 12:
          runp['Dfilename'] = dbase + '.dat.txt' + fend
          runp['Rfilename'] = rbase + '.ran.txt' + fend
        if DRopt == 13:
          runp['Dfilename'] = dbase + '.ran.txt' + fend
          runp['Rfilename'] = rbase + '.dat.txt' + fend
        if DRopt == 14:
          runp['Dfilename'] = dbase + '.ran.txt' + fend
          runp['Rfilename'] = rbase + '.ran.txt' + fend
         

      pfname = runp['foutbase'] + '-' + str(DRopt) + '.params'
      writexiparamfile(pfname,runp,Dftype=Dftype,Rftype=Rftype,DRopt=DRopt)
      mycmd = './xi %s' % (pfname)
      cmdlist.append(mycmd)
      if runopt == 0:
        ## write commands to std out, maybe later to a file for submission to a queue system?
        print mycmd 
      if runopt == 1:
        print 'running ',mycmd
        os.system(mycmd)
      if runopt == 2 and DRopt == 3:
        mycmd = 'nohup '+mycmd + ' &'
        print 'running ',mycmd
        os.system(mycmd)
      if runopt == 3 and DRopt < 3:
        print 'running ',mycmd
        os.system(mycmd)

  return cmdlist


  
## not sure this is useful -- see all the work in maskutils.py which is part of mksample
class mask:
  def __init__(self,fitsfname):
    self.mask = fitsio.read(fitsfname)

  def addradec(self,ax,color):
    myra = self.mask['RAMID'] + 90.0
    xx = np.where(myra > 360.0)
    myra[xx] -= 360.0

    myra = myra - 90.0

    ii = np.where(self.mask['WEIGHT'] > 0.01)[0]
    ax.plot(myra[ii],self.mask['DECMID'][ii],'.',color=color)

  def radecplot(self,prange=None,color='k'):
 
    f=plt.figure()
    ax=f.add_subplot(111)
    self.addradec(ax,color)
    if prange is not None:
      ax.axis(prange)

    return f, ax


