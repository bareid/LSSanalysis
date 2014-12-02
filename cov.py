import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm  ## for pcolormesh
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
  print 'this is not up to date!  i am using bootcov.py for this job'
  print 'BETH, you should edit/delete/merge for consistency and reduce defns of same functions in two places!!'
  return 0

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

def gettheorycov(taglist, tagname = 'v0', nsplits = 64, fbase = '/home/howdiedoo/boss/bethalexie/theoryerrors/makeHJcatalogv2zspace2', xitag='xiell',nboot=5000000):
  """
  taglist = labels of catalogs, since in practice it's nonsequential.  Currently completed list is
  101 - 123 and 201 - 229.
  nsplits = number of subboxes.  For current experiment, nsplits = 64.
  For comparison, we allow x number of averaging types to derive the cov.
  For production, we want avgtype = 0
  avgtype = 0: for each subbox region, average the statistic over all instances in the taglist.
  Then derive the cov from the variation among those subboxes.
  avgtype = 1: The average variation of a given realization about the mean of the subbox.
  Assumed format for the xi file labels:
  fbase + '_%d.%s.%d' % (tag,xitag,splitval)
  xitag = xiell or xiellcut or wp
  In case we change the taglist by adding more realizations, we will keep track of versions by tagname.
  Results: 
  """  
  for ss in range(nsplits):
    for tt in range(len(taglist)):
      tag = taglist[tt]
      if ss == 0 and tt == 0:
        svec, xiin = np.loadtxt(fbase + '_%d.%s.%d' % (tag,xitag,ss),unpack = True)
        nbins = len(xiin)
        ## allocate averages.
        xiavgsplit = np.zeros([nsplits,nbins])
        xiavgtot = np.zeros(nbins)
        ## allocate full storage of all the splits and tags individually.
        x = np.zeros([nsplits,len(taglist),nbins])
      else:
        xiin = np.loadtxt(fbase + '_%d.%s.%d' % (tag,xitag,ss),unpack = True,usecols=[1])
      x[ss,tt,:] = xiin[:]
      xiavgsplit[ss,:] = xiavgsplit[ss,:] + xiin[:]
      xiavgtot[:] = xiavgtot[:] + xiin[:]
    xiavgsplit[ss,:] = xiavgsplit[ss,:]/float(len(taglist))
  xiavgtot[:] = xiavgtot[:]/float(len(taglist))/float(nsplits)
  

      
  cov1 = np.zeros([nbins,nbins])
  cov2 = np.zeros([nbins,nbins])
  for i in range(nbins):
    for j in range(nbins):
      for ss in range(nsplits):
        cov1[i,j] += (xiavgsplit[ss,i] - xiavgtot[i])*(xiavgsplit[ss,j] - xiavgtot[j])
        for tt in range(len(taglist)):
          cov2[i,j] += (x[ss,tt,i] - xiavgsplit[ss,i]) * (x[ss,tt,j] - xiavgsplit[ss,j])

      cov1[i,j] = cov1[i,j]/float(nsplits)
      cov2[i,j] = cov2[i,j]/float(nsplits)/float(len(taglist))

  outfname1 = fbase + '_%s' % (tagname) + '.cov1'
  outfname2 = fbase + '_%s' % (tagname) + '.cov2'
  outfname3 = fbase + '_%s' % (tagname) + '.covboot'

  printcov(np.array(cov1),outfname1)
  printcov(np.array(cov2),outfname2)

  Cguess = np.zeros([nbins,nbins],dtype='float128')

  for b in range(nboot):
    rr = np.random.random_integers(0,nsplits-1,nsplits)
    xitrial = (xiavgsplit[rr,:]).sum(axis=0)/float(nsplits)
    xvec = np.matrix([xitrial-xiavgtot])
    Cguess += (xvec.T*xvec)

  Cguess = Cguess/float(nboot-1)
  Cguess = Cguess * nsplits ## regular bootstrap error give you error on the whole survey given chunks of the whole survey; here we want error on the small box so we need to multiply by ratio of the volumes = nsplits

  #### now let's compute icov for all these.
  ## eqn 17 of 0608064:
  p = len(Cguess[:,0])
  unbiasicov = float(nsplits - p - 2)/float(nsplits-1)
  print 'unbiasicov corr:',unbiasicov
  Cguess = np.matrix(Cguess,dtype='float64')
  invCguess = Cguess.I*unbiasicov
  printcov(np.array(Cguess),outfname3)


  #return svec, cov1, cov2
  ## tmp!
  return svec, cov1, cov2, xiavgsplit, xiavgtot, Cguess, invCguess

def gettheorycovxiMwp(taglist, tagname = 'v0', nsplits = 64, fbase = '/home/howdiedoo/boss/bethalexie/theoryerrors/makeHJcatalogv2zspace2', xitag='xiellcut',wptag = 'wp', nboot=5000000):
  """
  taglist = labels of catalogs, since in practice it's nonsequential.  Currently completed list is
  101 - 123 and 201 - 229.
  nsplits = number of subboxes.  For current experiment, nsplits = 64.
  For comparison, we allow x number of averaging types to derive the cov.
  For production, we want avgtype = 0
  avgtype = 0: for each subbox region, average the statistic over all instances in the taglist.
  Then derive the cov from the variation among those subboxes.
  avgtype = 1: The average variation of a given realization about the mean of the subbox.
  Assumed format for the xi file labels:
  fbase + '_%d.%s.%d' % (tag,xitag,splitval)
  xitag = xiell or xiellcut
  In case we change the taglist by adding more realizations, we will keep track of versions by tagname.
  Results: 
  """  
  for ss in range(nsplits):
    for tt in range(len(taglist)):
      tag = taglist[tt]
      if ss == 0 and tt == 0:
        svec, xiin = np.loadtxt(fbase + '_%d.%s.%d' % (tag,xitag,ss),unpack = True)
        nbins = len(xiin)
        ## read in wp too.
        rsig, wpin = np.loadtxt(fbase + '_%d.%s.%d' % (tag,wptag,ss),unpack = True)


        ## allocate averages.
        xiavgsplit = np.zeros([nsplits,nbins])
        xiavgtot = np.zeros(nbins)
        ## allocate full storage of all the splits and tags individually.
        x = np.zeros([nsplits,len(taglist),nbins])
      else:
        xiin = np.loadtxt(fbase + '_%d.%s.%d' % (tag,xitag,ss),unpack = True,usecols=[1])
      x[ss,tt,:] = xiin[:]
      xiavgsplit[ss,:] = xiavgsplit[ss,:] + xiin[:]
      xiavgtot[:] = xiavgtot[:] + xiin[:]
    xiavgsplit[ss,:] = xiavgsplit[ss,:]/float(len(taglist))
  xiavgtot[:] = xiavgtot[:]/float(len(taglist))/float(nsplits)
  

      
  cov1 = np.zeros([nbins,nbins])
  cov2 = np.zeros([nbins,nbins])
  for i in range(nbins):
    for j in range(nbins):
      for ss in range(nsplits):
        cov1[i,j] += (xiavgsplit[ss,i] - xiavgtot[i])*(xiavgsplit[ss,j] - xiavgtot[j])
        for tt in range(len(taglist)):
          cov2[i,j] += (x[ss,tt,i] - xiavgsplit[ss,i]) * (x[ss,tt,j] - xiavgsplit[ss,j])

      cov1[i,j] = cov1[i,j]/float(nsplits)
      cov2[i,j] = cov2[i,j]/float(nsplits)/float(len(taglist))

  outfname1 = fbase + '_%s' % (tagname) + '.cov1'
  outfname2 = fbase + '_%s' % (tagname) + '.cov2'
  outfname3 = fbase + '_%s' % (tagname) + '.covboot'

  printcov(np.array(cov1),outfname1)
  printcov(np.array(cov2),outfname2)

  Cguess = np.zeros([nbins,nbins],dtype='float128')

  for b in range(nboot):
    rr = np.random.random_integers(0,nsplits-1,nsplits)
    xitrial = (xiavgsplit[rr,:]).sum(axis=0)/float(nsplits)
    xvec = np.matrix([xitrial-xiavgtot])
    Cguess += (xvec.T*xvec)

  Cguess = Cguess/float(nboot-1)
  Cguess = Cguess * nsplits ## regular bootstrap error give you error on the whole survey given chunks of the whole survey; here we want error on the small box so we need to multiply by ratio of the volumes = nsplits

  #### now let's compute icov for all these.
  ## eqn 17 of 0608064:
  p = len(Cguess[:,0])
  unbiasicov = float(nsplits - p - 2)/float(nsplits-1)
  print 'unbiasicov corr:',unbiasicov
  Cguess = np.matrix(Cguess,dtype='float64')
  invCguess = Cguess.I*unbiasicov
  printcov(np.array(Cguess),outfname3)


  #return svec, cov1, cov2
  ## tmp!
  return svec, cov1, cov2, xiavgsplit, xiavgtot, Cguess, invCguess

def covadddummyzeros(covin,outfname,zlist):
  """
  add zeros to the cov matrix in front of current locations in cov.
  """
  nnew = len(zlist)
  try:
    (nx, ny) = covin.shape
  except:
    'wrong format for covin, doing nothing!'
    return None
  if nx != ny:
    'covin is not square, doing nothing!'
    return None
    
  cov = np.zeros([nx+nnew, nx+nnew])
  inew = 0
  zi = 0
  for iold in range(nx):
    if iold == zlist[zi]:
      inew = inew + 1
      if zi < len(zlist) - 1:
        zi += 1

    jnew = 0
    zj = 0
    for jold in range(nx):
      if jold == zlist[zj]:
      ## give zero rows a diagonal element so inversion works.
      ## does this work??
        if zlist[zi] < nx:
          cov[jnew,jnew] = covin[zlist[zi], zlist[zi]]
        else:
          cov[jnew,jnew] = covin[nx,nx]
        jnew = jnew + 1
        if zj < len(zlist) - 1:
          zj += 1
      ## copy old one into new location.
#      print 'inserting ',iold,jold,'into ',inew,jnew
      cov[inew,jnew] = covin[iold,jold]
      jnew += 1
    inew += 1
  for zi in range(len(zlist)):
    if zlist[zi] < nx:
      cov[zlist[zi]+zi, zlist[zi]+zi] = covin[zlist[zi], zlist[zi]] 
    else:
      cov[zlist[zi]+zi, zlist[zi]+zi] = covin[nx,nx]
  printcov(np.array(cov),outfname)
  return cov

def covremovedummyzeros(covin,zlist):
  """
  zlist is same as you submitted to add rows.
  """
  ## number of zeros that were added.
  nnew = len(zlist)
  assert (zlist[:-1] < zlist[1:]).all()
  try:
    (nx, ny) = covin.shape
  except:
    'wrong format for covin, doing nothing!'
    return None
  if nx != ny:
    'covin is not square, doing nothing!'
    return None
  x = copy.deepcopy(covin)
  for zz in zlist:
    x = np.delete(np.delete(x,zz,axis=0),zz,axis=1)

  return x

def icovsubsample(covin,zlist):
  """
  Read in original big covariance (not icov!).  
  zlist is the rows/cols you want to keep. 
  This routine subsamples, invert, fill back in with zeros all 
  along the INVERSE cov that's returned.
  """

  lbig = len(covin[:,0])
  ll = len(zlist)
  covnew = np.zeros([ll,ll])
  icovfinal = np.zeros([lbig, lbig])
  for i in range(ll):
    for j in range(ll):
      covnew[i,j] = covin[zlist[i], zlist[j]]

  icovnew = np.array(((np.matrix(covnew)).I))
  for i in range(ll):
    for j in range(ll):
      icovfinal[zlist[i], zlist[j]] = icovnew[i,j]

  return icovfinal

def covfillzeros(covin,keeplist):
  """
  Fill everything with zeros except for keeplist, then put 
  """

def plotdiagcov(ax,cov,scalevec=None,xvec=None,subindx=None,plotstr='k-'):
  """
  multiply the diagonal errors by scalevec, if you want to show errors on s*xi(s) for instance instead of xi.
  use subindx if you want to show only xi0 or only xi2 rather than their combination present in the cov.
  Plot diagerrs vs xvec if not None, otherwise just an index.
  """
  try:
    (nx, ny) = cov.shape
  except:
    'wrong format for covin, doing nothing!'
    return None
  if nx != ny:
    'covin is not square, doing nothing!'
    return None

  diagerr = np.array([(cov[i,i])**0.5 for i in range(nx)])
  if scalevec is not None:
    if len(scalevec) == len(diagerr):
      diagerr = diagerr*scalevec

  if subindx is not None:
    diagerr = diagerr[subindx]
    if scalevec is not None:
      if len(scalevec) == len(diagerr):
        diagerr = diagerr*scalevec

  if xvec is None:
    xvec = np.arange(0,len(diagerr)-0.5,1)
  ax.plot(xvec,diagerr,plotstr)
  
def getnormcov(cov,subindx1=None,subindx2=None):
  """
  plot the covariance matrix in subindx1 vs subindx2, normalized by the diagonal elements.
  normally called the correlation matrix.
  """
  try:
    (nx, ny) = cov.shape
  except:
    'wrong format for covin, doing nothing!'
    return None
  if nx != ny:
    'covin is not square, doing nothing!'
    return None
  if subindx1 is None:
    subindx1 = np.arange(0,nx-0.5,1,dtype='int')
  if subindx2 is None:
    subindx2 = np.arange(0,ny-0.5,1,dtype='int')

  nx2 = len(subindx1)
  ny2 = len(subindx2)
  corr = np.zeros([nx2,ny2])

  for i in range(nx2):
    for j in range(ny2):
      corr[i,j] = cov[subindx1[i],subindx2[j]]/np.sqrt(cov[subindx1[i],subindx1[i]]*cov[subindx2[j],subindx2[j]])
  
  return corr

def plotnormcov(cov,subindx1=None,subindx2=None,plotstr='k-',addcolorbar=True,ax=None,fsizex=4.34,fsizey=3.34):
  corr = getnormcov(cov,subindx1=subindx1,subindx2=subindx2)
  cmap = plt.get_cmap()
  levels=np.arange(-1.00,1.001,0.1)
  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
  if addcolorbar == True:
    assert ax is None
    plt.figure(figsize=[fsizex,fsizey])
  if ax is None:
    im = plt.pcolormesh(corr,cmap=cmap,norm=norm)
    if addcolorbar == True:
      plt.colorbar()
  else:
    ax.pcolormesh(corr,cmap=cmap,norm=norm)
  
  ### this doesn't work.  i think colorbar creates another axis, that's why this is a problem.
  #ax.pcolormesh(corr,cmap=cmap,norm=norm)
  #plt.colorbar()

def plotnormxi02cov(cov,nxi0=10,nxi2=10):
  """
  Need to make the figure sized by hand to get the color bar on there, for figure out how to make a separate ax for colorbar.
  """
  corr = getnormcov(cov,subindx1=None,subindx2=None)
  fig = plt.figure(figsize=[17,5])
  ax = fig.add_subplot(131,aspect='equal')
  plotnormcov(corr,subindx1=np.arange(0,nxi0-0.5,1,dtype='int'),subindx2=np.arange(0,nxi0-0.5,1,dtype='int'),addcolorbar=False,ax=ax)
  ax = fig.add_subplot(132,aspect='equal')
  plotnormcov(corr,subindx1=np.arange(nxi0,nxi0+nxi2-0.5,1,dtype='int'),subindx2=np.arange(nxi0,nxi0+nxi2-0.5,1,dtype='int'),addcolorbar=False,ax=ax)
  ax = fig.add_subplot(133,aspect='equal')
  plotnormcov(corr,subindx1=np.arange(0,nxi0-0.5,1,dtype='int'),subindx2=np.arange(nxi0,nxi0+nxi2-0.5,dtype='int'),addcolorbar=True,ax=None)
  return fig

def samplecovcorr(nbins,nparams,nsamples):
  A = 2./float((nsamples-nbins-1)*(nsamples-nbins-4))
  B = (nsamples-nbins-2.)/float((nsamples-nbins-1.)*(nsamples-nbins-4.))
  return (1.+B*(nbins-nparams))/float(1+A+B*(nparams+1))


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



