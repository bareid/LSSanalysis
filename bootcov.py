import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import re
import ximisc
import xiell
import wp
import xiwp
import os

## note I determined best splitxi0 and splitxi2 in 
#/Users/bareid/work/montserratdata/boss/tiledmockboss5002redo/comparetiledmockstotruthv0.ipynb
## for bin1.txt, it's what I put for defaults below: splitxi0=5,splitxi2=6)
## for bin1fineMU.txt (with rperpcut), it's splitxi0=1,splitxi2=1)



####### begin testing crap.
def tmpcompare():
  for i in range(200):
    f1 = "testing/testo%d" % i
    f2 = "/home/howdiedoo/boss/zdistvXlogbinsompcleverLSsmallscale/outputmksamplelatestdr10v7/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN-%03d_rmax48deltalog10rrebin-bin1.xielltrueNEW" % i
    print 'diff',i
    mystr = 'diff %s %s' % (f1, f2) 
    os.system(mystr)

def tmpcompare2():
  for i in range(200):
    f1 = "testing/testo%d" % i
    f2 = "/home/howdiedoo/boss/zdistvXlogbinsompcleverLSsmallscale/outputmksamplelatestdr10v7/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN-%03d_rmax48deltalog10rrebin-bin1fineMU.xiellcut" % i
    print 'diff',i
    mystr = 'diff %s %s' % (f1, f2) 
    os.system(mystr)

def tmpcompare3():
  for i in range(200):
    f1 = "testing/testo%d" % i
    f2 = "/home/howdiedoo/boss/zdistvXlogbinsompcleverLSsmallscale/outputmksamplelatestdr10v7/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN-%03d_rmax48deltalog10rrebin-bin1.xielltrueapproxRR" % i
    print 'diff',i
    mystr = 'diff %s %s' % (f1, f2) 
    os.system(mystr)
####### end testing crap.


## begin real useful stuff.

def checksymmetrycov(cov):
  for i in range(len(cov[:,0])):
    for j in range(len(cov[0,:])):
      assert cov[i,j] == cov[j,i]

def printcov(cov,fname):
  nx, ny = cov.shape
  ofp = open(fname,'w')
  for i in range(nx):
    for j in range(ny):
      ofp.write('%.12e ' % (cov[i,j]))
    ofp.write('\n')
  ofp.close()

def xicorrect(xiNNin, xiangin,splitxi0=5,splitxi2=6):
  xicorrxi=copy.deepcopy(xiangin.xi)
  xicorrxi[0,splitxi0:] = xiNNin.xi[0,splitxi0:]
  xicorrxi[1,splitxi2:] = xiNNin.xi[1,splitxi2:]
  xicorr = xiell.xiell(sxilist=[xiangin.svec,xicorrxi])
  ## need to fix xi0, xi2, xilong, xi. do those go through?
  return xicorr

def wpcorrect(wpNNin, wpangin, splitwp, wpstart, wpend):
  wpcorrwp=copy.deepcopy(wpangin.wp[wpstart:wpend+1])
  rsigin = copy.deepcopy(wpangin.rsig[wpstart:wpend+1])
  wpcorrwp[splitwp-wpstart:wpend+1-wpstart] = wpNNin.wp[splitwp:wpend+1]
  wpcorr = wp.wp(rpwplist=[rsigin,wpcorrwp])
  return wpcorr

def xiwpcorrect(xiNNin, xiangin,splitxi0,splitxi2,wpNNin, wpangin, splitwp, wpstart, wpend):
  """
  This function is for new statistic xiwp combining xi and wp.
  """
  mywp = wpcorrect(wpNNin, wpangin, splitwp, wpstart, wpend)
  myxiell = xicorrect(xiNNin, xiangin,splitxi0,splitxi2)
  myxiwp = xiwp.xiwp(myxiell,mywp)  
  return myxiwp

def debiasdataandcovwp(wpNNd, wpangd, wpangdhigh, wpangdlow, wpNNm, wpangm, wp012m, splitwp, wpstart,wpend,covstatfname,fname=None):
  """
  subtract the bias measured from the tiled mocks from the data, return a debiased combination.
  print it to a file (fname) to be fed to bethalexie code in long format.
  Also take in statistical covariance matrix and add two sources of systematics.
  """
  wpcorrdtmp = wpcorrect(wpNNd, wpangd, splitwp,wpstart,wpend)
  wpcorrm = wpcorrect(wpNNm, wpangm, splitwp,wpstart,wpend)
  wpdebiased = copy.deepcopy(wpcorrdtmp.wp)
  mydelta = wp012m.wp[wpstart:] - wpcorrm.wp
  print 'fractional wp correction:'
  print mydelta/wpcorrdtmp.wp
  wpdebiased = wpdebiased + mydelta

  ## now the cov.
  ## make sure this is the cov for the corrected statistic with same splits.
  if(0==0):
#  try:
    cov = np.loadtxt(covstatfname)
    assert len(cov[:,0]) == len(wpdebiased)
    splitz = covstatfname.split('splitswp')[1].split('_')
    assert len(splitz) >= 2
    ilist=[]
    for ss in splitz[:2]:
      ilist.append(int(ss))
    assert ilist[0] == splitwp
    assert ilist[1] == wpstart

    ## new jan 2 2014!!!  forgot to take into account the unbiasicov fac.  derive if from 
    ## product of cov and icov.
    ## guess icovfname
    tmp = covstatfname.split('/')
    tmp[-1] = 'i'+tmp[-1]
    icovstatfname = '/'.join(tmp)
    icov = np.loadtxt(icovstatfname)

    unbiasicovfac = (ximisc.getmatrixdiag(np.matrix(cov)*np.matrix(icov))).mean()
    print 'using htis unbiasicovfac correction, dividing cov by this',unbiasicovfac
    cov = cov/unbiasicovfac

    ndatacorr = len(wpcorrdtmp.wp)

    diagstat = np.zeros(ndatacorr)
    diagtot = np.zeros(ndatacorr)
    diagvar = np.zeros(ndatacorr)
    for i in range(len(diagstat)):
      diagstat[i] = cov[i,i]
    ## this must agree with wpcorrect assignmeents!
    wpangdiffvar = (0.5*(wpangdhigh.wp-wpangdlow.wp))**2 
    diagvar[0:splitwp-wpstart] = wpangdiffvar[wpstart:splitwp]
    print 'wp ang high/low variance: ',diagvar/diagstat
    print 'bias correction: ',mydelta
    diagvar = diagvar + (mydelta.flatten())**2
    print 'bias variance contribution: ',(mydelta.flatten())**2/diagstat
    ## add it into the covarianace matrix.
    for i in range(ndatacorr):
      cov[i,i] += diagvar[i]
      diagtot[i] = cov[i,i]
    print 'total sys variance fraction',diagtot/diagstat
  
    ## make it a matrix.
    cov = np.matrix(cov)
    icov = cov.I
  
    fcovout = covstatfname+'.sys'
    ## print the covariance and icov to new file.
    printcov(cov,fcovout)
    tmp = fcovout.split('/')
    tmp[-1] = 'i'+tmp[-1]
    ifcovout = '/'.join(tmp)
    printcov(icov,ifcovout)
  
    wpfinal = wp.wp(rpwplist=[wpcorrdtmp.rsig,wpdebiased],icovfname=ifcovout)    
    if fname is not None:
      wpfinal.printwp(fname)
    return wpfinal, cov

  else:
#  except:
    print 'cov file name does not match input splits, returning None!'
    wpfinal = wp.wp(rpwplist=[wpcorrdtmp.rsig,wpdebiased])    
    if fname is not None:
      wpfinal.printwp(fname)
    return wpfinal, None
  

## subtract the bias measured from the tiled mocks from the data, return a debiased combination.
## print it to a file to be fed to bethalexie code in long format.
def debiasdataandcov(xiNNd, xiangd, xiangdhigh, xiangdlow, xiNNm, xiangm, xi012m,splitxi0, splitxi2,covstatfname,nell=2,fname=None):
  """
  subtract the bias measured from the tiled mocks from the data, return a debiased combination.
  print it to a file (fname) to be fed to bethalexie code in long format.
  Also take in statistical covariance matrix and add two sources of systematics.
  """
  xicorrdtmp = xicorrect(xiNNd, xiangd, splitxi0,splitxi2)
  xicorrm = xicorrect(xiNNm, xiangm, splitxi0, splitxi2)
  xidebiased = copy.deepcopy(xicorrdtmp.xi)
  mydelta = xi012m.xi - xicorrm.xi
  xidebiased = xidebiased + mydelta

  ## now the cov.
  ## make sure this is the cov for the corrected statistic with same splits.
  if(0==0):
#  try:
    cov = np.loadtxt(covstatfname)
    assert len(cov[:,0]) == xiNNd.ndata
    splitz = covstatfname.split('splits')[1].split('_')
    assert len(splitz) == nell
    ilist=[]
    for ss in splitz:
      ilist.append(int(ss))
    assert ilist[0] == splitxi0
    assert ilist[1] == splitxi2

    ## new jan 2 2014!!!  forgot to take into account the unbiasicov fac.  derive if from 
    ## product of cov and icov.
    ## guess icovfname
    tmp = covstatfname.split('/')
    tmp[-1] = 'i'+tmp[-1]
    icovstatfname = '/'.join(tmp)
    icov = np.loadtxt(icovstatfname)

    unbiasicovfac = (ximisc.getmatrixdiag(np.matrix(cov)*np.matrix(icov))).mean()
    print 'using this unbiasicovfac correction, dividing cov by this',unbiasicovfac
    cov = cov/unbiasicovfac

    diagstat = np.zeros(xiNNd.ndata)
    diagtot = np.zeros(xiNNd.ndata)
    for i in range(len(diagstat)):
      diagstat[i] = cov[i,i]
    diagvar = np.zeros(xiNNd.ndata)
    xiangdiffvar = (0.5*(xiangdhigh.xi.flatten()-xiangdlow.xi.flatten()))**2 
    diagvar[0:splitxi0] = xiangdiffvar[0:splitxi0]
    nxi0 = len(xiNNd.xi0)
    diagvar[nxi0:nxi0+splitxi2] = xiangdiffvar[nxi0:nxi0+splitxi2] 
    print 'ang high/low variance: ',diagvar/diagstat
    diagvar = diagvar + (mydelta.flatten())**2
    print 'bias variance: ',(mydelta.flatten())**2/diagstat
    ## add it into the covarianace matrix.
    for i in range(xiNNd.ndata):
      cov[i,i] += diagvar[i]
      diagtot[i] = cov[i,i]
    print 'total sys variance fraction',diagtot/diagstat
  
    ## make it a matrix.
    cov = np.matrix(cov)
    icov = cov.I
  
    fcovout = covstatfname+'.sys'
    ## print the covariance and icov to new file.
    printcov(cov,fcovout)
    tmp = fcovout.split('/')
    tmp[-1] = 'i'+tmp[-1]
    ifcovout = '/'.join(tmp)
    printcov(icov,ifcovout)
  
    xifinal = xiell.xiell(sxilist=[xiNNd.svec,xidebiased],icovfname=ifcovout)
    if fname is not None:
      ofp = open(fname,'w')
      ofp.write("# ellmax = %d\n" % ((nell-1)*2))
      for i in range(len(xifinal.svec.flatten())):
        ofp.write('%e %e\n' % (xifinal.svec.flatten()[i], xifinal.xi.flatten()[i]))
      ofp.close()

    return xifinal, cov

  else:
#  except:
    print 'cov file name does not match input splits, returning None!'
    xifinal = xiell.xiell(sxilist=[xiNNd.svec,xidebiased])
    if fname is not None:
      ofp = open(fname,'w')
      ofp.write("# ellmax = %d\n" % ((nell-1)*2))
      for i in range(len(xifinal.svec.flatten())):
        ofp.write('%e %e\n' % (xifinal.svec.flatten()[i], xifinal.xi.flatten()[i]))
      ofp.close()
      return xifinal, None
    
def debiasdataandcovxiMwp(xiNNd, xiangd, xiangdhigh, xiangdlow, xiNNm, xiangm, xi012m,splitxi0, splitxi2, wpNNd, wpangd, wpangdhigh, wpangdlow, wpNNm, wpangm, wp012m, splitwp, wpstart,wpend,covstatfname,nell=2,fname=None):
#def debiasdataandcovwp(wpNNd, wpangd, wpangdhigh, wpangdlow, wpNNm, wpangm, wp012m, splitwp, wpstart,wpend,covstatfname,fname=None):
  """
  subtract the bias measured from the tiled mocks from the data, return a debiased combination.
  print it to a file (fname) to be fed to bethalexie code in long format.
  Also take in statistical covariance matrix and add two sources of systematics.
  """
#def xiwpcorrect(xiNNin, xiangin,splitxi0,splitxi2,wpNNin, wpangin, splitwp, wpstart, wpend):
  xiwpcorrdtmp = xiwpcorrect(xiNNd, xiangd, splitxi0,splitxi2,\
                         wpNNd, wpangd, splitwp, wpstart, wpend)

  xiwpcorrm = xiwpcorrect(xiNNm, xiangm, splitxi0, splitxi2,\
                        wpNNm,wpangm,splitwp,wpstart,wpend)

  xiwpdebiased = copy.deepcopy(xiwpcorrdtmp)
  #tmp!
#  print xiwpdebiased.xiell
#  print xiwpdebiased.wp

  mydeltaxi = xi012m.xi - xiwpcorrm.xiell.xi  ## subtract xi objects.
  mydeltawp = wp012m.wp[wpstart:wpend+1] - xiwpcorrm.wp.wp
  xiwpdebiased.xiell.xi = xiwpdebiased.xiell.xi + mydeltaxi
  xiwpdebiased.wp.wp = xiwpdebiased.wp.wp + mydeltawp
  xiwpdebiased.xiwp = np.concatenate((xiwpdebiased.xiell.xilong, xiwpdebiased.wp.wp))
  xiwpanghigh = xiwp.xiwp(xiangdhigh,wpangdhigh)
  xiwpanglow = xiwp.xiwp(xiangdlow,wpangdlow)

  ## now the cov.
  ## make sure this is the cov for the corrected statistic with same splits.
  if(0==0):
#  try:
    cov = np.loadtxt(covstatfname)
    assert len(cov[:,0]) == xiwpdebiased.ntot
    splitz = covstatfname.split('splitswp')[0].split('splits')[1].split('_')
    assert len(splitz) >= nell
    ilist=[]
    tmp=0
    for ss in splitz:
      tmp += 1
      ilist.append(int(ss))
      if tmp >= nell:
        break
    assert ilist[0] == splitxi0
    assert ilist[1] == splitxi2
    splitzwp = covstatfname.split('splitswp')[1].split('_')
    print splitzwp
    assert len(splitzwp) >= 3
    ilist=[]
    tmp=0
    for ss in splitzwp:
      tmp += 1
      ilist.append(int(ss))
      if tmp >= 3:
        break
    assert ilist[0] == splitwp
    assert ilist[1] == wpstart
    assert ilist[2] == wpend

    ## new jan 2 2014!!!  forgot to take into account the unbiasicov fac.  derive if from 
    ## product of cov and icov.
    ## guess icovfname
    tmp = covstatfname.split('/')
    tmp[-1] = 'i'+tmp[-1]
    icovstatfname = '/'.join(tmp)
    icov = np.loadtxt(icovstatfname)

    unbiasicovfac = (ximisc.getmatrixdiag(np.matrix(cov)*np.matrix(icov))).mean()
    print 'using this unbiasicovfac correction, dividing cov by this',unbiasicovfac
    cov = cov/unbiasicovfac

    diagstat = np.zeros(xiwpdebiased.ntot)
    diagtot = np.zeros(xiwpdebiased.ntot)
    diagvar = np.zeros(xiwpdebiased.ntot)
    for i in range(len(diagstat)):
      diagstat[i] = cov[i,i]

    xiangdiffvar = (0.5*(xiangdhigh.xi.flatten()-xiangdlow.xi.flatten()))**2 
    diagvar[0:splitxi0] = xiangdiffvar[0:splitxi0]
    nxi0 = len(xiNNd.xi0)
    nxi2 = len(xiNNd.xi2)
    diagvar[nxi0:nxi0+splitxi2] = xiangdiffvar[nxi0:nxi0+splitxi2] 
    wpangdiffvar = (0.5*(wpangdhigh.wp-wpangdlow.wp))**2
    diagvar[nxi0+nxi2:nxi0+nxi2+splitwp-wpstart] = wpangdiffvar[wpstart:splitwp]
    print 'ang high/low variance: ',diagvar/diagstat
    diagvar[0:nxi0+nxi2] = diagvar[0:nxi0+nxi2] + (mydeltaxi.flatten())**2
    diagvar[nxi0+nxi2:] = diagvar[nxi0+nxi2:] + (mydeltawp)**2
    print 'bias variance xi: ',(mydeltaxi.flatten())**2/diagstat[0:nxi0+nxi2]
    print 'bias variance wp: ',(mydeltawp)**2/diagstat[nxi0+nxi2:]
    ## add it into the covarianace matrix.
    for i in range(xiwpdebiased.ntot):
      cov[i,i] += diagvar[i]
      diagtot[i] = cov[i,i]
    print 'total sys variance fraction',diagtot/diagstat
  
    ## make it a matrix.
    cov = np.matrix(cov)
    icov = cov.I
  
    fcovout = covstatfname+'.sys'
    ## print the covariance and icov to new file.
    printcov(cov,fcovout)
    tmp = fcovout.split('/')
    tmp[-1] = 'i'+tmp[-1]
    ifcovout = '/'.join(tmp)
    printcov(icov,ifcovout)

    xiwpfinal = xiwp.xiwp(xiwpdebiased.xiell, xiwpdebiased.wp, icovfname=ifcovout)
    if fname is not None:
      ofp = open(fname,'w')
      ofp.write("# ellmax = %d\n" % ((nell-1)*2))
      for i in range(len(xiwpfinal.xiell.svec.flatten())):
        ofp.write('%e %e\n' % (xiwpfinal.xiell.svec.flatten()[i], xiwpfinal.xiell.xi.flatten()[i]))
      for i in range(len(xiwpfinal.wp.wp)):
        ofp.write('%e %e\n' % (xiwpfinal.wp.rsig[i], xiwpfinal.wp.wp[i]))
      ofp.close()

    return xiwpfinal, cov

  else:
    print 'cov file name does not match input splits, returning None!'
    xiwpfinal = xiwp.xiwp(xiwpdebiased.xiell, xiwpdebiased.wp, icovfname=ifcovout)
    if fname is not None:
      ofp = open(fname,'w')
      ofp.write("# ellmax = %d\n" % ((nell-1)*2))
      for i in range(len(xiwpfinal.xi.svec.flatten())):
        ofp.write('%e %e\n' % (xiwpfinal.xiell.svec.flatten()[i], xiwpfinal.xiell.xi.flatten()[i]))
      for i in range(len(xiwpfinal.wp.wp)):
        ofp.write('%e %e\n' % (xiwpfinal.wp.rsig[i], xiwpfinal.wp.wp[i]))
      ofp.close()

    return xiwpfinal, None

def parsebootinfo(bootfile,workingdir):
  ## stuff we need to get from the file.
  nsub = None
  pixelfname = None
  #nsubdir = None
  fbaseNNstart = None
  fbaseangstart = None
  ## end stuff.

  ifp = open(bootfile,'r')
  for line in ifp:
    if(re.match('nsub:',line)):
      nsub = int(line.split('nsub:')[1].strip('\n').strip(' '))
    if(re.match('pixelfname:',line)):
      pixelfname = workingdir + line.split('pixelfname:')[1].strip('\n').strip(' ')
#    if(re.match('nsubdir:',line)):
#      nsubdir = line.split('nsubdir:')[1].strip('\n').strip(' ')
    if(re.match('fbaseNNstart:',line)):
      fbaseNNstart = workingdir + line.split('fbaseNNstart:')[1].strip('\n').strip(' ')
    if(re.match('fbaseangstart:',line)):
      fbaseangstart = workingdir + line.split('fbaseangstart:')[1].strip('\n').strip(' ')


    if(re.match('fbaseNNtotN:',line)):
      fbaseNNtotN = workingdir + line.split('fbaseNNtotN:')[1].strip('\n').strip(' ')
    if(re.match('fbaseNNtotS:',line)):
      fbaseNNtotS = workingdir + line.split('fbaseNNtotS:')[1].strip('\n').strip(' ')
    if(re.match('fbaseangtotN:',line)):
      fbaseangtotN = workingdir + line.split('fbaseangtotN:')[1].strip('\n').strip(' ')
    if(re.match('fbaseangtotS:',line)):
      fbaseangtotS = workingdir + line.split('fbaseangtotS:')[1].strip('\n').strip(' ')
#  print nsub, pixelfname, fbaseNNstart, fbaseangstart, fbaseNNtotN, fbaseNNtotS, fbaseangtotN, fbaseangtotS
  return nsub, pixelfname, fbaseNNstart, fbaseangstart, fbaseNNtotN, fbaseNNtotS, fbaseangtotN, fbaseangtotS
      
def getpixlist(pixelfname,nsub):
  ## read in pixels.
  Pin = np.loadtxt(pixelfname)
  assert len(Pin) == nsub
  ## define dtype we need here.
  pixlist = np.zeros(nsub,dtype=[('PID','int'),('idec','int'),('ramin','float'),('ramax','float'),('decmin','float'),('decmax','float'),('NorS','int')])
  pixlist['PID'] = Pin[:,0]
  pixlist['NorS'] = Pin[:,1]
  pixlist['idec'] = Pin[:,2]
  pixlist['ramin'] = Pin[:,4]
  pixlist['ramax'] = Pin[:,5]
  pixlist['decmin'] = Pin[:,6]
  pixlist['decmax'] = Pin[:,7]
  return pixlist

## need to get DRfac and fixRRdown from N and S, out of the files, send to xiellfromDR.
## copy /home/howdiedoo/boss/bootstrapdr10v7/calcxi02bootcov.py for how to deal with some linear combination of NN and
## ang when deriving cov.
def getbootcov(bootfile, workingdir, covtag=None, NSortot=2, nboot = 5000000, fbaseend='_rmax48deltalog10r',\
               xiellorwp=0,rpimax=80.,splitwp=7,wpstart=1,wpend=19,\
               nell=3,binfile=None,rperpcut=-1.,smallRRcut=-1.,\
               dfacs=1,dfacmu=1,icovfname=None,smincut=-1.,smaxcut=1.e12,\
               splitxi0=5,splitxi2=6,fbaseendxiell='_rmax48deltalog10r',fbaseendwp='_xigrid'):
  """
  Get covariance matrix.
  fbaseend = '_rmax48deltalog10r' for xiell or '_xigrid' for wp.
  Third tier of stuff goes directly to xiellfromDR
  expect splitxi0/splitxi2 [those go to xicorrect; values determined in
  comparetiledcmockstotruthv0
  Added functionality for wp: xiellorwp = 1, splitwp = where to go from ang to NN.
  rpimax is for wp, default is 80.
  Leaving variable fbaseend for backward compatibility,  but if xiellorwp == 2, defaults to 
  using fbaseendxiell and fbaseendwp
  """
  #NNorang = 0 [NN] or 1 [ang] or 2 [optimal unbiased combination, not yet written]
  ## nevermind, that doesnt get used anywhere!!?? deleted, was the 4th elt in the list.

  assert xiellorwp >= 0 and xiellorwp <= 2

  nsub, pixelfname, fbaseNNstart, fbaseangstart, \
  fbaseNNtotN, fbaseNNtotS, fbaseangtotN, fbaseangtotS =  parsebootinfo(bootfile,workingdir)

  if nsub is None or pixelfname is None or fbaseNNstart is None or fbaseangstart is None:
    print 'bad boot file, getbootcov returning None!'
    return None
  pixlist = getpixlist(pixelfname,nsub)

  myfbase_NN = fbaseNNstart
  myfbase_ang = fbaseangstart

  if xiellorwp == 0 or xiellorwp == 1:
    DRfacN_NN, fixRRdownN_NN = ximisc.getDRfactors(fbaseNNtotN+fbaseend)
    DRfacS_NN, fixRRdownS_NN = ximisc.getDRfactors(fbaseNNtotS+fbaseend)
    DRfacN_ang, fixRRdownN_ang = ximisc.getDRfactors(fbaseangtotN+fbaseend)
    DRfacS_ang, fixRRdownS_ang = ximisc.getDRfactors(fbaseangtotS+fbaseend)
  else: ##xiwp statistic.  xiellorwp == 2
    ## xiell
    DRfacN_NNxiell, fixRRdownN_NNxiell = ximisc.getDRfactors(fbaseNNtotN+fbaseendxiell)
    DRfacS_NNxiell, fixRRdownS_NNxiell = ximisc.getDRfactors(fbaseNNtotS+fbaseendxiell)
    DRfacN_angxiell, fixRRdownN_angxiell = ximisc.getDRfactors(fbaseangtotN+fbaseendxiell)
    DRfacS_angxiell, fixRRdownS_angxiell = ximisc.getDRfactors(fbaseangtotS+fbaseendxiell)
    ## wp
    DRfacN_NNwp, fixRRdownN_NNwp = ximisc.getDRfactors(fbaseNNtotN+fbaseendwp)
    DRfacS_NNwp, fixRRdownS_NNwp = ximisc.getDRfactors(fbaseNNtotS+fbaseendwp)
    DRfacN_angwp, fixRRdownN_angwp = ximisc.getDRfactors(fbaseangtotN+fbaseendwp)
    DRfacS_angwp, fixRRdownS_angwp = ximisc.getDRfactors(fbaseangtotS+fbaseendwp)


  if xiellorwp == 0:
    splittag = 'splits%d_%d' % (splitxi0,splitxi2)
  if xiellorwp == 1:
    splittag = 'splitswp%d_%d_%d' % (splitwp,wpstart,wpend)
  if xiellorwp == 2:
    splittagxiell = 'splits%d_%d' % (splitxi0,splitxi2)
    splittagwp = 'splitswp%d_%d_%d' % (splitwp,wpstart,wpend)
    splittag = splittagxiell+'_' + splittagwp

  if binfile is not None:
    bintag = binfile.split('/')[-1].split('.')[0]
    covoutNN = 'covtotv7NN_b%d_N%d_rebin-%s' % (nboot,nsub,bintag)
    covoutang = 'covtotv7ang_b%d_N%d_rebin-%s' % (nboot,nsub,bintag)
    covoutcorr = 'covtotv7corr_b%d_N%d_rebin-%s_%s' % (nboot,nsub,bintag,splittag)
  else:
    covoutNN = 'covtotv7NN_b%d_N%d' % (nboot,nsub)
    covoutang = 'covtotv7ang_b%d_N%d' % (nboot,nsub)
    covoutcorr = 'covtotv7corr_b%d_N%d_%s' % (nboot,nsub,splittag)
    

  if covtag is not None:
    covoutNN = covoutNN + '_%s' % covtag
    covoutang = covoutang + '_%s' % covtag
    covoutcorr = covoutcorr + '_%s' % covtag

  icovoutNN = 'i'+covoutNN
  icovoutang = 'i'+covoutang
  icovoutcorr = 'i'+covoutcorr

  
  if xiellorwp == 0 or xiellorwp == 1:
    DRinfoN_NN = [DRfacN_NN, fixRRdownN_NN]
    DRinfoS_NN = [DRfacS_NN, fixRRdownS_NN]
    DRinfoN_ang = [DRfacN_ang, fixRRdownN_ang]
    DRinfoS_ang = [DRfacS_ang, fixRRdownS_ang]
  else:
    #xiell
    DRinfoN_NNxiell = [DRfacN_NNxiell, fixRRdownN_NNxiell]
    DRinfoS_NNxiell = [DRfacS_NNxiell, fixRRdownS_NNxiell]
    DRinfoN_angxiell = [DRfacN_angxiell, fixRRdownN_angxiell]
    DRinfoS_angxiell = [DRfacS_angxiell, fixRRdownS_angxiell]
    # wp
    DRinfoN_NNwp = [DRfacN_NNwp, fixRRdownN_NNwp]
    DRinfoS_NNwp = [DRfacS_NNwp, fixRRdownS_NNwp]
    DRinfoN_angwp = [DRfacN_angwp, fixRRdownN_angwp]
    DRinfoS_angwp = [DRfacS_angwp, fixRRdownS_angwp]

  for ns in range(nsub):
    print ns
    fbase_NN = myfbase_NN + ('-%03d' % (ns))+fbaseend
    fbase_ang = myfbase_ang + ('-%03d' % (ns))+fbaseend
    if xiellorwp == 2:
      fbase_NNxiell = myfbase_NN + ('-%03d' % (ns))+fbaseendxiell
      fbase_angxiell = myfbase_ang + ('-%03d' % (ns))+fbaseendxiell
      fbase_NNwp = myfbase_NN + ('-%03d' % (ns))+fbaseendwp
      fbase_angwp = myfbase_ang + ('-%03d' % (ns))+fbaseendwp

    xx = np.where(pixlist['PID'] == ns)[0]
    assert len(xx) == 1
    assert xx[0] == ns
    NorSval = pixlist['NorS'][xx[0]]
    if(NorSval == 0):
      if xiellorwp == 0 or xiellorwp == 1:
        DRinfo_NN = DRinfoN_NN
        DRinfo_ang = DRinfoN_ang
      else:
        DRinfo_NNxiell = DRinfoN_NNxiell
        DRinfo_angxiell = DRinfoN_angxiell
        DRinfo_NNwp = DRinfoN_NNwp
        DRinfo_angwp = DRinfoN_angwp

    else: #south
      if xiellorwp == 0 or xiellorwp == 1:
        DRinfo_NN = DRinfoS_NN
        DRinfo_ang = DRinfoS_ang
      else:
        DRinfo_NNxiell = DRinfoS_NNxiell
        DRinfo_angxiell = DRinfoS_angxiell
        DRinfo_NNwp = DRinfoS_NNwp
        DRinfo_angwp = DRinfoS_angwp


    if xiellorwp == 0:
      xiinNN = xiell.xiellfromDR(fbase_NN,nell,binfile,rperpcut,dfacs,dfacmu,icovfname,smincut,smaxcut,DRinfo_NN,smallRRcut)
      xiinang = xiell.xiellfromDR(fbase_ang,nell,binfile,rperpcut,dfacs,dfacmu,icovfname,smincut,smaxcut,DRinfo_ang,smallRRcut)
      xicorr = xicorrect(xiinNN, xiinang, splitxi0, splitxi2)
    if xiellorwp == 1:  ## doing wp
      xiinNNtmp = wp.wpfromDR(fbase_NN,DRfacinfo=DRinfo_NN,rpimax=rpimax,icovfname=icovfname)
      xiinangtmp = wp.wpfromDR(fbase_ang,DRfacinfo=DRinfo_ang,rpimax=rpimax,icovfname=icovfname)

      ## these are for later, saving cov of NN and ang separately.
      xiinNN = wp.wpfromDR(fbase_NN,DRfacinfo=DRinfo_NN,rpimax=rpimax,icovfname=icovfname,wpstart=wpstart,wpend=wpend)
      xiinang = wp.wpfromDR(fbase_ang,DRfacinfo=DRinfo_ang,rpimax=rpimax,icovfname=icovfname,wpstart=wpstart,wpend=wpend)

      ##wpstart,end not already applied to this NN and ang!
      xicorr = wpcorrect(xiinNNtmp,xiinangtmp,splitwp,wpstart,wpend)

    if xiellorwp == 2:  ##doing xiwp
      xiinNNxiell = xiell.xiellfromDR(fbase_NNxiell,nell,binfile,rperpcut,dfacs,dfacmu,icovfname,smincut,smaxcut,DRinfo_NNxiell,smallRRcut)
      xiinangxiell = xiell.xiellfromDR(fbase_angxiell,nell,binfile,rperpcut,dfacs,dfacmu,icovfname,smincut,smaxcut,DRinfo_angxiell,smallRRcut)

      xiinNNwptmp = wp.wpfromDR(fbase_NNwp,DRfacinfo=DRinfo_NNwp,rpimax=rpimax,icovfname=icovfname)
      xiinangwptmp = wp.wpfromDR(fbase_angwp,DRfacinfo=DRinfo_angwp,rpimax=rpimax,icovfname=icovfname)

      xiinNNwp = wp.wpfromDR(fbase_NNwp,DRfacinfo=DRinfo_NNwp,rpimax=rpimax,icovfname=icovfname,wpstart=wpstart,wpend=wpend)
      xiinangwp = wp.wpfromDR(fbase_angwp,DRfacinfo=DRinfo_angwp,rpimax=rpimax,icovfname=icovfname,wpstart=wpstart,wpend=wpend)
      xiinNN = xiwp.xiwp(xiinNNxiell,xiinNNwp)
      xiinang = xiwp.xiwp(xiinangxiell,xiinangwp)
      xicorr = xiwpcorrect(xiinNNxiell, xiinangxiell,splitxi0,splitxi2,xiinNNwptmp, xiinangwptmp, splitwp, wpstart, wpend)

    ## tmp!  we tested to make sure we recovered the same correlation fxns as with old code.  Good!
    #tmpfname = "testing/testo%d" % ns
    #xiin.printxiellshort(tmpfname)
    if(ns == 0):
      if(xiellorwp == 0):
        ndata = xiinNN.ndata
        ndatacorr = ndata
      if(xiellorwp == 1):
        ndata = len(xiinNN.wp)
        ndatacorr = len(xicorr.wp)
      if(xiellorwp == 2):
        ndata = xiinNN.ntot
        ndatacorr = ndata

      xilistNN = np.zeros([nsub,ndata],dtype='float128')
      xilistang = np.zeros([nsub,ndata],dtype='float128')
      xilistcorr = np.zeros([nsub,ndatacorr],dtype='float128')
    if(xiellorwp == 0):
      xilistNN[ns,:] = xiinNN.xilong
      xilistang[ns,:] = xiinang.xilong
      xilistcorr[ns,:] = xicorr.xilong
    if(xiellorwp == 1):
      xilistNN[ns,:] = xiinNN.wp
      xilistang[ns,:] = xiinang.wp
      xilistcorr[ns,:] = xicorr.wp
    if(xiellorwp == 2):
      xilistNN[ns,:] = xiinNN.xiwp
      xilistang[ns,:] = xiinang.xiwp
      xilistcorr[ns,:] = xicorr.xiwp

  ## check means with total counts.
  nindx = np.where(pixlist['NorS'] == 0)[0]
  sindx = np.where(pixlist['NorS'] == 1)[0]
  print 'N/S: ',len(nindx), len(sindx)

  ## now compute mean and bootstrap errors:
  if(NSortot == 0):
    ximeanNN = (xilistNN[nindx,:]).sum(axis=0)/float(len(nindx))
    ximeanang = (xilistang[nindx,:]).sum(axis=0)/float(len(nindx))
    ximeancorr = (xilistcorr[nindx,:]).sum(axis=0)/float(len(nindx))
    ntot = len(nindx)
    ## restrict xilist to N only
    xilistNN = xlistNN[nindx,:]
    xilistang = xlistang[nindx,:]
    xilistcorr = xlistcorr[nindx,:]

  if(NSortot == 1):
    ximeanNN = (xilistNN[sindx,:]).sum(axis=0)/float(len(sindx))
    ximeanang = (xilistang[sindx,:]).sum(axis=0)/float(len(sindx))
    ximeancorr = (xilistcorr[sindx,:]).sum(axis=0)/float(len(sindx))
    ntot = len(sindx)
    ## restrict xilist to S only
    xilistNN = xlistNN[sindx,:]
    xilistang = xlistang[sindx,:]
    xilistcorr = xlistcorr[sindx,:]

  if(NSortot == 2):
    ximeanNN = xilistNN.sum(axis=0)/float(nsub)
    ximeanang = xilistang.sum(axis=0)/float(nsub)
    ximeancorr = xilistcorr.sum(axis=0)/float(nsub)
    ntot = nsub
 
  xitotNN = np.zeros(ndata,dtype='float128')
  xitotang = np.zeros(ndata,dtype='float128')
  xitotcorr = np.zeros(ndatacorr,dtype='float128')
  CguessNN = np.zeros([ndata,ndata],dtype='float128')
  Cguessang = np.zeros([ndata,ndata],dtype='float128')
  Cguesscorr = np.zeros([ndatacorr,ndatacorr],dtype='float128')

  for b in range(nboot):
    rr = np.random.random_integers(0,ntot-1,ntot)
    xitrialNN = (xilistNN[rr,:]).sum(axis=0)/float(ntot)
    xitrialang = (xilistang[rr,:]).sum(axis=0)/float(ntot)
    xitrialcorr = (xilistcorr[rr,:]).sum(axis=0)/float(ntot)
    xvecNN = np.matrix([xitrialNN-ximeanNN])
    xvecang = np.matrix([xitrialang-ximeanang])
    xveccorr = np.matrix([xitrialcorr-ximeancorr])
    CguessNN += (xvecNN.T*xvecNN)
    Cguessang += (xvecang.T*xvecang)
    Cguesscorr += (xveccorr.T*xveccorr)

  CguessNN = CguessNN/float(nboot-1)
  Cguessang = Cguessang/float(nboot-1)
  Cguesscorr = Cguesscorr/float(nboot-1)

  ## put this back in after tests.
  #### now let's compute icov for all these.
  ## eqn 17 of 0608064:
  p = len(CguessNN[:,0])
  unbiasicov = float(ntot - p - 2)/float(ntot-1)

  CguessNN = np.matrix(CguessNN,dtype='float64')
  invCguessNN = CguessNN.I*unbiasicov 
  printcov(CguessNN,covoutNN)
  printcov(invCguessNN,icovoutNN)

  Cguessang = np.matrix(Cguessang,dtype='float64')
  invCguessang = Cguessang.I*unbiasicov 
  printcov(Cguessang,covoutang)
  printcov(invCguessang,icovoutang)

  Cguesscorr = np.matrix(Cguesscorr,dtype='float64')
  invCguesscorr = Cguesscorr.I*unbiasicov 
  printcov(Cguesscorr,covoutcorr)
  printcov(invCguesscorr,icovoutcorr)

  return CguessNN, invCguessNN, Cguessang, invCguessang, Cguesscorr, invCguesscorr

def covaddsys(covfname,splitxi0=5,splitxi2=6):
  """
  input statistical "corr" covariance matrix filename.
  we infer splitxi0/2 from the filename.
  """
  ## check that the splits match the statistical cov file names.
  ## BETHHERE!! 

def getpixlistcolors(pixlist,clist):
  """
  Blah
  """
  mycsorted = {}
  pixlistcpy = copy.deepcopy(pixlist)
  ## don't screw up original list
  pixlistcpy.sort(order=('idec','ramin'))
  npix = len(pixlist['PID'])
  nn = 0
  ncolors = len(clist)
  icstart = 0
  ic = icstart
  mycsorted[pixlistcpy['PID'][nn]] = clist[ic % ncolors]
  curridec = pixlistcpy['idec'][nn]
  ic += 1
  nn += 1
  while nn < npix:
    while pixlistcpy['idec'][nn] == curridec:
      mycsorted[pixlistcpy['PID'][nn]] = clist[ic % ncolors]
      ic += 1
      nn += 1
      if nn == npix:
        break
    if nn == npix:
      break
    ## otherwise, new row!
    icstart = icstart + 1
    ic = icstart
    mycsorted[pixlistcpy['PID'][nn]] = clist[ic % ncolors]
    curridec = pixlistcpy['idec'][nn]
    ic += 1
    nn += 1
      
  return mycsorted ## returns dictionary of color for each pixel value.

## plot bootstrap regions.
## copying from ~/boss/bootstrap/plotbootregionsworphans.py
def plotbootstrapregions(bootfile, workingdir,ax=None,lw=3,plotstr='k-',clist=['b.','g.','r.','c.','m.'],Doutbase='/home/howdiedoo/boss/mksamplecatslatestdr10/v7/threed/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN.txt'):
  """
  Blah.
  """
  nclist = len(clist)

  if ax is None:
    ff = plt.figure(figsize=[6,6])
    ax=ff.add_subplot(1,1,1)
  else:
    ff = None


  nsub, pixelfname, fbaseNNstart, fbaseangstart, \
  fbaseNNtotN, fbaseNNtotS, fbaseangtotN, fbaseangtotS =  parsebootinfo(bootfile,workingdir)

  if nsub is None or pixelfname is None or fbaseNNstart is None or fbaseangstart is None:
    print 'bad boot file, getbootcov returning None!'
    return None
  pixlist = getpixlist(pixelfname,nsub)
  mycdict = getpixlistcolors(pixlist,clist)

  ## convert back to degrees.
  radtodeg = 180./np.pi
  pixlist['ramin'] = pixlist['ramin']*radtodeg
  pixlist['ramax'] = pixlist['ramax']*radtodeg
  pixlist['decmin'] = pixlist['decmin']*radtodeg
  pixlist['decmax'] = pixlist['decmax']*radtodeg



  for ns in range(nsub):
    ### bring back after I set up colors.
    ffD = Doutbase+'.'+str(ns)
    raD, decD = np.loadtxt(ffD,unpack=True,usecols=[0,1])
    raD = raD - 90.0
    xxD = np.where(raD < 0.)[0]
    if(len(xxD) > 0):
      raD[xxD] += 360.
    plt.plot(raD,decD,mycdict[ns])

  for ns in range(nsub):
    plt.plot([pixlist['ramin'][ns], pixlist['ramin'][ns]], [pixlist['decmin'][ns], pixlist['decmax'][ns]],plotstr,linewidth=lw)
    plt.plot([pixlist['ramax'][ns], pixlist['ramax'][ns]], [pixlist['decmin'][ns], pixlist['decmax'][ns]],plotstr,linewidth=lw)
    plt.plot([pixlist['ramin'][ns], pixlist['ramax'][ns]], [pixlist['decmin'][ns], pixlist['decmin'][ns]],plotstr,linewidth=lw)
    plt.plot([pixlist['ramin'][ns], pixlist['ramax'][ns]], [pixlist['decmax'][ns], pixlist['decmax'][ns]],plotstr,linewidth=lw)

  return ff, ax

def plotbootstrapcumhist(bootfile,workingdir,ax=None,fcatbase='/home/howdiedoo/boss/mksamplecatslatestdr10/v7/threed/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN'):

  Doutbase = fcatbase + '.txt'
  Routbase = fcatbase + '.ran.txt'

  if ax is None:
    ff = plt.figure(figsize=[6,6])
    ax=ff.add_subplot(1,1,1)
  else:
    ff = None

  nsub, pixelfname, fbaseNNstart, fbaseangstart, \
  fbaseNNtotN, fbaseNNtotS, fbaseangtotN, fbaseangtotS =  parsebootinfo(bootfile,workingdir)

  if nsub is None or pixelfname is None or fbaseNNstart is None or fbaseangstart is None:
    print 'bad boot file, getbootcov returning None!'
    return None
  pixlist = getpixlist(pixelfname,nsub)

  ## convert back to degrees.
  radtodeg = 180./np.pi
  pixlist['ramin'] = pixlist['ramin']*radtodeg
  pixlist['ramax'] = pixlist['ramax']*radtodeg
  pixlist['decmin'] = pixlist['decmin']*radtodeg
  pixlist['decmax'] = pixlist['decmax']*radtodeg

  nDlist = []
  nRlist = []

  for ns in range(nsub):
    ### bring back after I set up colors.
    ffD = Doutbase+'.'+str(ns)
    ffR = Routbase+'.'+str(ns)
    #raD, decD = np.loadtxt(ffD,unpack=True,usecols=[0,1])
    #raR, decR = N.loadtxt(ffR,unpack=True,usecols=[0,1])
    raD, wgtD = np.loadtxt(ffD,unpack=True,usecols=[0,3])
    raR, wgtR = np.loadtxt(ffR,unpack=True,usecols=[0,3])
    nDlist.append(float(wgtD.sum()))
    nRlist.append(float(wgtR.sum()))

  return np.array(nDlist), np.array(nRlist)
    
  



if __name__ == '__main__':
  #parsebootinfo('/home/howdiedoo/boss/bootstrapdr10v7/bootNsub200.dat','/home/howdiedoo/boss/')
  parsebootinfo('/home/howdiedoo/boss/bootstrapdr10v7/bootNsub200.dat','/home/howdiedoo/boss/')
  workingdir = '/home/howdiedoo/boss/'
  workingdir = '/home/howdiedoo/boss/'
  bootfile = workingdir + 'bootstrapdr10v7/bootNsub200.dat'
  binfile = workingdir + 'zdistvXlogbinsompcleverLSsmallscale/bin1.txt'
  binfile2 = workingdir + 'zdistvXlogbinsompcleverLSsmallscale/bin1fineMU.txt'
  covtag = 'testo'

  ##passed!
  #getbootcov(bootfile, workingdir, covtag, nboot = 50000, binfile=binfile)
  #print 'running tmpcompare now!'
  #tmpcompare()

  ##passed!
  #getbootcov(bootfile, workingdir, covtag, nboot = 50000, binfile=binfile,smallRRcut=400.)
  #print 'running tmpcompare3 now for approxRR!'
  #tmpcompare3()


  ## now we want to compare xiellcut.
  getbootcov(bootfile, workingdir, covtag, nboot = 50000, binfile=binfile2,rperpcut=5.336546e-01)
  print 'running tmpcompare2 now for xiellcut!'
  tmpcompare2()
  ## passed!  there are tiny differences still in inner bin (most are much smaller)
  ## anyway, they are much different than the statistical errors.  But i will print the bin file
  ## just in case.
#< 7.852356e-01 1.632387e+01 -2.219600e+00 -1.323648e+01
#> 7.852356e-01 1.643968e+01 -1.921517e+00 -1.340287e+01



