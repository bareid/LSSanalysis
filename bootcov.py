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
import boot

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

def printmean(dmean,fname):
  ofp = open(fname,'w')
  for i in range(len(dmean)):
    ofp.write('%.12e\n' % (dmean[i]))
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

def xiwpvec(xiellin,wpin,wpcrossin,wpstart,wpend):
      #dvec = xiwpvec(xitmp,wptmp,wpcrosstmp,wpstart,wpend)
  assert wpend+1 <= len(wpcrossin.wp) ## want to check before I do the below.
  myxi = np.concatenate((xiellin.xi0[1:], xiellin.xi2[1:])) ## skip first elt of xi0 and xi2 by default.
  if wpend+1 > len(wpcrossin.wp):
    return np.concatenate((myxi,wpcrossin.wp[wpstart:],wp.wp[len(wpcrossin):wpend+1]))
  else:
    return np.concatenate((myxi,wpcrossin.wp[wpstart:wpend+1]))


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

def parsebootinfo(bootfile):
  """
  Assumes current structure of mksamplecatslatestdr12 (input as basedir)
  separate output directories for all the different statistics.
  within the output directories, the subcat statistics are in nsubdir
  Rewritten for DR12.
  """
  ## stuff we need to get from the file.
  nsub = None
  nsubdir = None
  pixelfname = None
  fbase = None
  fbasetotN = None
  fbasetotS = None
  ## end stuff.

  ifp = open(bootfile,'r')
  for line in ifp:
    if(re.match('nsub:',line)):
      nsub = int(line.split('nsub:')[1].strip('\n').strip(' '))
    if(re.match('pixelfname:',line)):
      pixelfname = line.split('pixelfname:')[1].strip('\n').strip(' ')
    if(re.match('nsubdir:',line)):
      nsubdir = line.split('nsubdir:')[1].strip('\n').strip(' ')
#    if(re.match('nsubdir:',line)):
#      nsubdir = line.split('nsubdir:')[1].strip('\n').strip(' ')
    if(re.match('fbase:',line)):
      fbase = line.split('fbase:')[1].strip('\n').strip(' ')
    if(re.match('fbasetotN:',line)):
      fbasetotN = line.split('fbasetotN:')[1].strip('\n').strip(' ')
    if(re.match('fbasetotS:',line)):
      fbasetotS = line.split('fbasetotS:')[1].strip('\n').strip(' ')
    

  return nsub, nsubdir, pixelfname, fbase, fbasetotN, fbasetotS


## need to get DRfac and fixRRdown from N and S, out of the files, send to xiellfromDR.
## copy /home/howdiedoo/boss/bootstrapdr10v7/calcxi02bootcov.py for how to deal with some linear combination of NN and
## ang when deriving cov.
def getbootcov(bootfile, basedir, outdirbase = 'outputdr12', covoutfname=None, NSortot=2, nboot = 5000000, \
               rpimax=80.,wpstart=1,wpend=19,\
               nell=3,rperpcut=-1.,smallRRcut=-1.,\
               dfacs=1,dfacmu=1,icovfname=None,smincut=-1.,smaxcut=1.e12,\
               binfname_xiell='xibinfiles/bin1fineMU.txt',\
               nbar2d=[-1.,-1.],nbar3d=[-1.,-1],\
               whichtask=4):
## resurrect these later.
#               splitxi0=5,splitxi2=6,splitwp=7):
  """
  Get covariance matrix.
  We're going to do all tasks at once by default (4).
  whichtask = 0: xiell
  whichtask = 1: wp (compute xi(rp,rpi))
  whichtask = 2: wtheta
  whichtask = 3: Hogg spec-im cross-correlation.
  whichtask = 4: combine xiell and wp in usual way.

  Third tier of stuff goes directly to xiellfromDR
  rpimax is for wp, default is 80.
  nbar2d,nbar3d needs to be computed separately for N and S.
  """

  nsub, nsubdir, pixelfname, fbase, fbasetotN, fbasetotS = parsebootinfo(bootfile=basedir+bootfile)

  NSlist = [0,1]
  NStaglist = ['N','S']

  for xx in [nsub, nsubdir, pixelfname, fbase, fbasetotN, fbasetotS]:
    if xx is None:
      print 'bad bootfile!'
      return None

  b = boot.bootpix()
  b.readregions(basedir + pixelfname)
  assert b.nsub == nsub

  ## this list will be filled
  DRinfolist = [-1,-1,-1,-1]

  taglist= ['-xiell','-xigrid','-wtheta','-wpcross']

  ## get global DR factors for taglist.
  for ii in range(len(taglist)-1):
    tag = taglist[ii]
    tmp = np.zeros([2,2]) # first index is N or S.  DRfac, fixRR stored for each.

    for NS, NStag, ff in zip(NSlist, NStaglist,[fbasetotN,fbasetotS]):
      try:
      #if 0==0:
        tmp[NS,0], tmp[NS,1] = ximisc.getDRfactors(basedir + '/'+outdirbase + tag +'/'+ff)
      except:
        tmp[NS,:] = -1.
    DRinfolist[ii] = tmp.copy()

  ## now get DR info for wpcross.
  ### nevermind! we reduce this to two ratios.
  ## DRinfolist[3] = np.zeros([2,4,2])
  DRinfolist[3] = np.zeros([2,2])
  tag = taglist[3]
  for NS, NStag, ff in zip(NSlist, NStaglist,[fbasetotN,fbasetotS]):
    try:
      normfac = ximisc.getDRnormswpcross(basedir + '/'+outdirbase + tag +'/'+ff) 
      DRinfolist[3][NS][0] = normfac[0,0]/normfac[2,0]
      DRinfolist[3][NS][1] = normfac[0,1]/normfac[1,1]
    except:
      DRinfolist[3][NS][:] = -1.

  tasklist = np.zeros(4,dtype='int')
  if whichtask == 4:
    tasklist = np.array([1,1,0,1],dtype='int')
  else:
    tasklist[whichtask] = 1

  if tasklist[3] > 0:
    assert (nbar2d[:] > 0).all()
    assert (nbar3d[:] > 0).all()
    assert (DRinfolist[3][:,:].flatten() > 0).all()

  for ns in range(nsub):
    xx = np.where(b.pixlist['PID'] == ns)[0]
    assert len(xx) == 1
    assert xx[0] == ns
    NorSval = b.pixlist['NorS'][xx[0]]
    for tt in range(len(tasklist)):
      if tasklist[tt] == 0: continue
      tag = taglist[tt]
      ff = basedir+'/'+outdirbase + tag +'/' + nsubdir + '/' + fbase + '.%04d.Np' % (ns) 
      if tt == 0: #xiell
        xitmp = xiell.xiellfromDR(ff,binfile=binfname_xiell,rperpcut=rperpcut,nell=nell,smallRRcut=smallRRcut,dfacs=dfacs,dfacmu=dfacmu,smincut=smincut,smaxcut=smaxcut,DRfacinfo=DRinfolist[tt][NorSval]) 
        dvec = xitmp.xilong
      if tt == 1: #wp
        wptmp = wp.wpfromDR(ff,DRfacinfo=DRinfolist[tt][NorSval],rpimax=rpimax)
        dvec = wptmp.wp

      if tt == 2: #wtheta
        wttmp = wtheta.wthetafromDR(ff,DRfacinfo=DRinfolist[tt][NorSval])
        dvec = wttmp.wtheta

      if tt == 3: #wpcross
        wpcrosstmp = wp.wpcrossHogg(ff,DRfacinfo=DRinfolist[tt][NorSval],nbar2d=nbar2d[NorSval],nbar3d=nbar3d[NorSval])
        dvec = wpcrosstmp.wp

    if whichtask == 4:
      dvec = xiwpvec(xitmp,wptmp,wpcrosstmp,wpstart,wpend)

    if ns == 0:  ## allocate!
      ndata = len(dvec)
      dveclist = np.zeros([nsub,ndata],dtype='float128')
    dveclist[ns,:] = dvec[:]

  ## check means with total counts.
  nindx = np.where(b.pixlist['NorS'] == 0)[0]
  sindx = np.where(b.pixlist['NorS'] == 1)[0]
  nsindx = np.where((b.pixlist['NorS'] == 0) | (b.pixlist['NorS'] == 1))[0]
  print 'N/S: ',len(nindx), len(sindx), len(nsindx)
  assert len(nsindx) == nsub
  assert (nsindx == np.arange(0,nsub,1,dtype='int')).all()
  myindx= nsindx
  ## assume we want nsindx for this, but can restore N/S option later if I want.

  dmean = (dveclist[myindx,:]).sum(axis=0)/float(len(myindx))
  ntot = len(myindx)
  ntotflt = float(ntot)

  print 'hi beth'
  print dmean

  Cmat = np.zeros([ndata,ndata],dtype='float128')
  for b in range(nboot):
    rr = np.random.random_integers(0,ntot-1,ntot)
    dtrial = (dveclist[rr,:]).sum(axis=0)/ntotflt
    xvec = np.matrix([dtrial-dmean])
    Cmat += (xvec.T*xvec)

  Cmat = Cmat/float(nboot-1)
  Cmat = np.matrix(Cmat,dtype='float64')
  iCmat = Cmat.I ##
  print 'not assuming any bootstrap unbias factor for now!'
  if covoutfname is not None:
    printcov(Cmat,covoutfname)
    printcov(iCmat,covoutfname+'.inv')
    printmean(dmean,covoutfname+'.mean')
  return Cmat, iCmat, dmean









