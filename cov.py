import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import re
import ximisc
import xiell
import os

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

def checksymmetrycov(cov):
  """
  if you're sending a matrix object, bracket with np.array()
  """
  for i in range(len(cov[:,0])):
    for j in range(len(cov[0,:])):
      assert cov[i,j] == cov[j,i]

def printcov(cov,fname):
  """
  if you're sending a matrix object, bracket with np.array()
  """
  ofp = open(fname,'w')
  for i in range(len(cov[:,0])):
    for j in range(len(cov[0,:])):
      ofp.write('%.12e ' % (cov[i,j]))
    ofp.write('\n')
  ofp.close()

## need to get DRfac and fixRRdown from N and S, out of the files, send to xiellfromDR.
## copy /home/howdiedoo/boss/bootstrapdr10v7/calcxi02bootcov.py for how to deal with some linear combination of NN and
## ang when deriving cov.
def getbootcov(bootfile, workingdir, covtag, NNorang=0, NSortot=2, nboot = 5000000, fbaseend='_rmax48deltalog10r',\
               nell=3,binfile=None,rperpcut=-1.,smallRRcut=-1.,
               dfacs=1,dfacmu=1,icovfname=None,smincut=-1.,smaxcut=1.e12):
  """
  Get covariance matrix.
  fbaseend = '_rmax48deltalog10r' for xiell or '_xigrid' for wp.
  NNorang = 0 [NN] or 1 [ang] or 2 [optimal unbiased combination, not yet written]
  Third tier of stuff goes directly to xiellfromDR
  """
  nsub, pixelfname, fbaseNNstart, fbaseangstart, \
  fbaseNNtotN, fbaseNNtotS, fbaseangtotN, fbaseangtotS =  parsebootinfo(bootfile,workingdir)

  if nsub is None or pixelfname is None or fbaseNNstart is None or fbaseangstart is None:
    print 'bad boot file, getbootcov returning None!'
    return None
  pixlist = getpixlist(pixelfname,nsub)

  if NNorang == 0:
    myfbase = fbaseNNstart
    DRfacN, fixRRdownN = ximisc.getDRfactors(fbaseNNtotN+fbaseend)
    DRfacS, fixRRdownS = ximisc.getDRfactors(fbaseNNtotS+fbaseend)

  elif NNorang == 1:
    myfbase = fbaseangstart
    DRfacN, fixRRdownN = ximisc.getDRfactors(fbaseangtotN+fbaseend)
    DRfacS, fixRRdownS = ximisc.getDRfactors(fbaseangtotS+fbaseend)
  else:
    print 'NNorang = ',NNorang,'not supported.'
    return None

  DRinfoN = [DRfacN, fixRRdownN]
  DRinfoS = [DRfacS, fixRRdownS]

  for ns in range(nsub):
    fbase = myfbase + ('-%03d' % (ns))+fbaseend
    xx = np.where(pixlist['PID'] == ns)[0]
    assert len(xx) == 1
    assert xx[0] == ns
    NorSval = pixlist['NorS'][xx[0]]
    if(NorSval == 0):
      DRinfo = DRinfoN
    else:
      DRinfo = DRinfoS
    xiin = xiell.xiellfromDR(fbase,nell,binfile,rperpcut,dfacs,dfacmu,icovfname,smincut,smaxcut,DRinfo,smallRRcut)
    ## tmp!  we tested to make sure we recovered the same correlation fxns as with old code.  Good!
    tmpfname = "testing/testo%d" % ns
    xiin.printxiellshort(tmpfname)
    if(ns == 0):
      ndata = xiin.ndata
      xilist = np.zeros([nsub,ndata],dtype='float128')
    xilist[ns,:] = xiin.xilong

  ## check means with total counts.
  nindx = np.where(pixlist['NorS'] == 0)[0]
  sindx = np.where(pixlist['NorS'] == 1)[0]
  print 'N/S: ',len(nindx), len(sindx)

  ## now compute mean and bootstrap errors:
  if(NSortot == 0):
    ximean = (xilist[nindx,:]).sum(axis=0)/float(len(nindx))
    ntot = len(nindx)
    ## restrict xilist to N only
    xilist = xlist[nindx,:]
  if(NSortot == 1):
    ximean = (xilist[sindx,:]).sum(axis=0)/float(len(sindx))
    ntot = len(sindx)
    xilist = xlist[sindx,:]
  if(NSortot == 2):
    ximean = xilist.sum(axis=0)/float(nsub)
    ntot = nsub

  xitot = np.zeros(ndata,dtype='float128')
  Cguess = np.zeros([ndata,ndata],dtype='float128')

  for b in range(nboot):
    rr = np.random.random_integers(0,ntot-1,ntot)
    xitrial = (xilist[rr,:]).sum(axis=0)/float(ntot)
    xvec = np.matrix([xitrial-ximean])
    Cguess += (xvec.T*xvec)

  Cguess = Cguess/float(nboot-1)

  #### now let's compute icov for all these.
  ## eqn 17 of 0608064:
  p = len(Cguess[:,0])
  unbiasicov = float(ntot - p - 2)/float(ntot-1)
  Cguess = np.matrix(Cguess,dtype='float64')
  invCguess = Cguess.I*unbiasicov 

#  printcov(Cguess,"cov.tmp")
#  printcov(invCguess,"icov.tmp")
  return Cguess, invCguess

def invcovsubsample(icov,onedi):
  """
  inverts icov to get cov, subsamples elts of cov to the ones you want, and reinverts.  
  Passes back cov, icov pair.
  onedi should be a numpy array of the elements of cov you want to keep.  They must be ordered/sorted and unique!
  Unless you're mixing up your data vector order compared to the cov.  We will not handle that case for now.
  """
  icov = np.matrix(icov)
  tmpcov = np.array(icov.I)
  ## make sure the elements of onedi are ordered.
  assert (onedi[:-1] < onedi[1:]).all()
  assert len(tmpcov[:,0]) >= len(onedi)
  mx, my = np.meshgrid(onedi,onedi)
  newcov = np.matrix(tmpcov[(mx,my)])
  newicov = newcov.I
  assert newcov.shape == (len(onedi), len(onedi))
  assert newicov.shape == (len(onedi), len(onedi))
  return newcov, newicov

def tmpcompare():
  for i in range(200):
    f1 = "testing/testo%d" % i
    f2 = "/Users/bareid/work/montserratdata/boss/zdistvXlogbinsompcleverLSsmallscale/outputmksamplelatestdr10v7/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN-%03d_rmax48deltalog10rrebin-bin1.xielltrueNEW" % i
    print 'diff',i
    mystr = 'diff %s %s' % (f1, f2) 
    os.system(mystr)

def tmpcompare2():
  for i in range(200):
    f1 = "testing/testo%d" % i
    f2 = "/Users/bareid/work/montserratdata/boss/zdistvXlogbinsompcleverLSsmallscale/outputmksamplelatestdr10v7/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN-%03d_rmax48deltalog10rrebin-bin1.xiellcut" % i
    print 'diff',i
    mystr = 'diff %s %s' % (f1, f2) 
    os.system(mystr)

def tmpcompare3():
  for i in range(200):
    f1 = "testing/testo%d" % i
    f2 = "/Users/bareid/work/montserratdata/boss/zdistvXlogbinsompcleverLSsmallscale/outputmksamplelatestdr10v7/dr10v7bootworphansNsub200/collidedBR-collate-cmass-dr10v7-FBBRNN-%03d_rmax48deltalog10rrebin-bin1.xielltrueapproxRR" % i
    print 'diff',i
    mystr = 'diff %s %s' % (f1, f2) 
    os.system(mystr)

if __name__ == '__main__':
  #parsebootinfo('/Users/bareid/work/montserratdata/boss/bootstrapdr10v7/bootNsub200.dat','/home/howdiedoo/boss/')
  parsebootinfo('/Users/bareid/work/montserratdata/boss/bootstrapdr10v7/bootNsub200.dat','/home/howdiedoo/boss/')
  workingdir = '/home/howdiedoo/boss/'
  workingdir = '/Users/bareid/work/montserratdata/boss/'
  bootfile = workingdir + 'bootstrapdr10v7/bootNsub200.dat'
  binfile = workingdir + 'zdistvXlogbinsompcleverLSsmallscale/bin1.txt'
  binfile2 = workingdir + 'zdistvXlogbinsompcleverLSsmallscale/bin1fineMU.txt'
  covtag = 'testo'

  getbootcov(bootfile, workingdir, covtag, nboot = 50000, binfile=binfile)
  print 'running tmpcompare now!'
  tmpcompare()

  getbootcov(bootfile, workingdir, covtag, nboot = 50000, binfile=binfile,smallRRcut=400.)
  print 'running tmpcompare3 now for approxRR!'
  tmpcompare3()


  ## now we want to compare xiellcut.
  getbootcov(bootfile, workingdir, covtag, nboot = 50000, binfile=binfile2,rperpcut=5.336546e-01)
  print 'running tmpcompare2 now for xiellcut!'
  tmpcompare2()



