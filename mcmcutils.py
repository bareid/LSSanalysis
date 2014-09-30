import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
#import os
import re
import hod
import BRutils

def wgtresult(wgt,v):
  wgtsum = wgt.sum()
  m = (wgt*v).sum()/wgtsum
  s = ((wgt*v**2).sum()/wgtsum - m**2)**0.5
  return m, s

def wgtresult2(wgt,v,flist=[0.682689492137,0.954499736104]):

  vcpy = v.copy()
  wgtcpy = wgt.copy()
  xx = vcpy.argsort()
  vcpy = vcpy[xx]
  wgtcpy = wgtcpy[xx]
  wgtsum = wgtcpy.sum()
  xxi = np.where(wgtcpy.cumsum() > wgtsum*0.5)[0][0]
  vcen = vcpy[xxi]
  result = [vcen]
  for ff in flist:
    vlow = vcpy[np.where(wgtcpy.cumsum() > wgtsum*(0.5-ff*0.5))[0][0]]
    vhigh = vcpy[np.where(wgtcpy.cumsum() > wgtsum*(0.5+ff*0.5))[0][0]]
    result.append(vcen-vlow)
    result.append(vhigh-vcen)
  return result

def wgtresult2oneside(wgt,v,flist=[0.682689492137,0.954499736104],loworhigh=0):

  vcpy = v.copy()
  wgtcpy = wgt.copy()
  xx = vcpy.argsort()
  vcpy = vcpy[xx]
  wgtcpy = wgtcpy[xx]
  wgtsum = wgtcpy.sum()
  result = []
  for ff in flist:
    if loworhigh == 0:
      v1d = vcpy[np.where(wgtcpy.cumsum() > wgtsum*ff)[0][0]]
    else:
      v1d = vcpy[np.where(wgtcpy.cumsum() > wgtsum*(1.-ff))[0][0]]
    result.append(v1d)
  return result

def wgtresult2raw(wgt,v,flist=[0.682689492137,0.954499736104]):

  vcpy = v.copy()
  wgtcpy = wgt.copy()
  xx = vcpy.argsort()
  vcpy = vcpy[xx]
  wgtcpy = wgtcpy[xx]
  wgtsum = wgtcpy.sum()
  xxi = np.where(wgtcpy.cumsum() > wgtsum*0.5)[0][0]
  vcen = vcpy[xxi]
  result = [vcen]
  for ff in flist:
    vlow = vcpy[np.where(wgtcpy.cumsum() > wgtsum*(0.5-ff*0.5))[0][0]]
    vhigh = vcpy[np.where(wgtcpy.cumsum() > wgtsum*(0.5+ff*0.5))[0][0]]
    result.append(vlow)
    result.append(vhigh)
  return np.array(result)

def wgtresultint(wgt,v,flist=[0.682689492137,0.954499736104]):
  result = wgtresult2(wgt,v,flist)
  return result[0], 0.5*(result[1]+result[2])

def getclevels(Hin,flist=[0.683,0.954]):
  """
  Input a list of fractions, output the histogram values
  where we should cut them.
  """
  clist = []
  d = Hin.copy()
  ## -1 means "whatever is needed; only giving one number means i want a 1d array.
  d.shape = -1
  d.sort()
  
  for ff in flist:
# = d.reshape(len(myhflip[0,:])*len(myhflip([:,0]),
    target = 1.-ff
    xx = np.where(d.cumsum()/float(d.sum()) < target)[0]
    yy = np.where(d.cumsum()/float(d.sum()) > target)[0]

    clev = 0.5*(d[xx[-1]]+d[yy[0]])
#    print d.cumsum()[xx[-1]], d.cumsum()[yy[0]], target
    xx = np.where(Hin > clev)
    clist.append(clev)
  return clist

def contourplot(x,y,wgt=None,ax=None,\
                 nxbins=25,nybins=25,flist=[0.683,0.954],
                 linestyle='-',color='k',linewidths=3):
  """
  Input optional weight to make a weighted 2d histogram.
  flist is fraction of points you want to enclose.
  """
  if wgt is None:
    H, xedges, yedges = np.histogram2d(x,y,[nxbins,nybins])
  else:
    H, xedges, yedges = np.histogram2d(x,y,[nxbins,nybins],weights=wgt)
  Hflip = np.zeros([nybins,nxbins])  # for some idiotic reason, contour progam bins are flipped.
  xcen = (xedges[1:] + xedges[:-1])*0.5
  ycen = (yedges[1:] + yedges[:-1])*0.5
  for x in range(nxbins):
    for y in range(nybins):
      Hflip[y,x] = H[x,y]

  clist = getclevels(Hflip,flist)
  if ax is None:
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.contour(xcen,ycen,Hflip,levels=clist,linestyles=linestyle,colors=color,linewidths=linewidths)
    return f, ax
  else:
    ax.contour(xcen,ycen,Hflip,levels=clist,linestyles=linestyle,colors=color,linewidths=linewidths)    

## can delete this?? just duplicated with contourplot above.  oops!
def contourplotgeneric(vx,vy,wgt=None,ax=None,wgtopt=1,\
               nxbins=25,nybins=25,flist=[0.683,0.954],
               linestyle='-',color='k',linewidths=3):
  """
  Set wgtopt = 1 to weight the points according to the field 'weight'
  in the chain.  Set wgtopt = 0 to weight the points equally.
  flist is fraction of chain you want to enclose.
  """
  if wgt is None:
    wgtopt = 0
  if wgtopt == 1:  
    H, xedges, yedges = np.histogram2d(vx,vy,[nxbins,nybins],weights=wgt)
  else:
    H, xedges, yedges = np.histogram2d(vx,vy,[nxbins,nybins])

  Hflip = np.zeros([nybins,nxbins])  # for some idiotic reason, contour progam bins are flipped.
  xcen = (xedges[1:] + xedges[:-1])*0.5
  ycen = (yedges[1:] + yedges[:-1])*0.5
  for x in range(nxbins):
    for y in range(nybins):
      Hflip[y,x] = H[x,y]

  clist = getclevels(Hflip,flist)
  if ax is None:
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.contour(xcen,ycen,Hflip,levels=clist,linestyles=linestyle,colors=color,linewidths=linewidths)
    return f, ax
  else:
    ax.contour(xcen,ycen,Hflip,levels=clist,linestyles=linestyle,colors=color,linewidths=linewidths)


class chain:
  def __init__(self,chainfname,colfname,takelog10m=0):
    """
    To homogenize chains, take the log10 of mass columns if takelog10m = 1
    This function also specifies mcmc parameters and which are apparently varying.
    """

    names = []
    ifpp = open(colfname,'r')
    mcmcp = {} ## dictionary relating column names to index in a mcmc parameter list.
    mcmcpreverse = {} ## dictionary going from integer index to name.
    mcmcpcnt = 0

    for line in ifpp:
      for nn in line.split(','):
        mynn = nn.strip(' ').strip('\n')
        names.append(mynn)
        if re.match('weight',mynn) or re.match('chi2',mynn) or \
           re.match('lnlike',mynn) or \
           re.match('nbar',mynn) or re.match('fsat',mynn):
          continue  ## not an mcmc parameter.
        mcmcp[mynn] = mcmcpcnt
        mcmcpreverse[mcmcpcnt] = mynn
        mcmcpcnt += 1
    self.mcmcp = mcmcp
    self.mcmcpreverse = mcmcpreverse

    try:
      cc = np.genfromtxt(chainfname,names=names)
    except:
      print 'misamtch between %s and %s.  Returning None!' % (chainfname,names)
    if takelog10m == 1:
      cc['M_min'][:] = np.log10(cc['M_min'][:])
      cc['M1'][:] = np.log10(cc['M1'][:])
      cc['M_cut'][:] = np.log10(cc['M_cut'][:])

    self.chain = cc
    self.names = names
    (self.nelts, self.ncols) = (len(cc)), len(names)
    self.fname = chainfname

    ## figure out which parameters are fixed and which are varying.
    ## this appears to work!
    mcmcfixed = np.zeros(len(mcmcp),dtype='int')
    for i in range(len(mcmcpreverse)):
      m, std = (cc[mcmcpreverse[i]][:]).mean(), (cc[mcmcpreverse[i]][:]).std()
      if(std == 0.):
        mcmcfixed[i] = 1
      ## deal with rounding errors
      if(np.fabs(m) > 2.0e-6):
        if std/np.fabs(m) < 2.0e-6:
          mcmcfixed[i] = 1
      else:
        if std < 2.0e-6:
          mcmcfixed[i] = 1
    self.mcmcfixed = mcmcfixed


  def __str__(self):
    mystr = 'chain of shape: '+str(self.chain.shape)+'\n'
    mystr += 'chain imported from %s\n' % self.fname
    mystr += 'chain elts: %s\n' % self.names
    return mystr

  def minchi2(self):
    mymin = (self.chain['chi2_tot'][:]).min()
    xx = np.where(self.chain['chi2_tot'][:] == mymin)
    return xx, mymin



  ## copying from plot2dhistcleanboth2planelv5talk.py in prettyplotsv2
  def contourplot(self,xname,yname,ax=None,wgtopt=1,\
                 nxbins=25,nybins=25,flist=[0.683,0.954],
                 linestyle='-',color='k',linewidths=3):
    """
    Set wgtopt = 1 to weight the points according to the field 'weight'
    in the chain.  Set wgtopt = 0 to weight the points equally.
    flist is fraction of chain you want to enclose.
    """
    wgt = None
    if wgtopt == 1:
      wgt = self.chain['weight']
    x = self.chain[xname]
    y = self.chain[yname]
    contourplot(x,y,wgt,ax,nxbins,nybins,flist,linestyle,color,linewidths)


  
  def fillstepmatrix(self,steprescale=2.4):
    """
    Generate step matrix from the mcmc chain, noting which parameters are actually being varied.
    Multiply step_mat by steprescale/np.sqrt(npvary)
    std steprescale is 2.4.
    """

    nparam = len(self.mcmcp)
    meanlist = np.zeros(nparam)
    wgtsum = (self.chain['weight'][:]).sum()

    ## calculate the mean for every parameter.
    for p in range(nparam):
      nn = self.mcmcpreverse[p]  ## name of parameter p
      if(self.mcmcfixed[p] == 1):
        meanlist[p] = self.chain[nn][0]
        #assert (self.chain[nn][:]).std() == 0.
      else:
        meanlist[p] = (self.chain['weight'][:] * self.chain[nn][:]).sum()/wgtsum
    ## calculate covariance matrix for the parameters that are varying.
    yy = np.where(self.mcmcfixed == 0)[0]  ## yy contains indices of varying parameters.
    npvary = len(yy)
    cov = np.zeros([npvary, npvary])
    for p in range(npvary):
      pi = yy[p]
      pinn = self.mcmcpreverse[yy[p]]  ## name index in original list of mcmc parameters, including non-varying parameters.
      for q in range(npvary):
        qi = yy[q]
        qinn = self.mcmcpreverse[yy[q]]
        cov[p,q] = (self.chain[pinn][:] * self.chain[qinn][:] * self.chain['weight'][:]).sum()/wgtsum - \
                   meanlist[pi] * meanlist[qi]

    covmat = np.matrix(cov)
    #generate step matrix.  copying from MW's file cov.py in mcmc_wprp
    eigval, eigvec = np.linalg.eig(covmat)
    eigval = np.real(eigval)
    ww   = np.nonzero(eigval<1e-10)[0]

    if len(ww)>0:
      print "Found some worrying eigenvalues:"
      print eigval
      print "remapping those <1e-10 to zero."
      eigval[ww] = 0
      print eigval

    step_mat = np.zeros([npvary,npvary])
    for p in range(npvary):
      for q in range(npvary):
        for a in range(npvary):
          step_mat[p,q] += eigvec[p,a]*np.sqrt(eigval[a])*eigvec[q,a]

    # Scale by an optimization factor and write the result.
    step_mat *= steprescale/np.sqrt((npvary)*1.0)
    step_mat_big = np.zeros([nparam, nparam])
    

    for p in range(npvary):
      pi = yy[p]
      pinn = self.mcmcpreverse[yy[p]]  ## name index in original list of mcmc parameters, including non-varying parameters.
      for q in range(npvary):
        qi = yy[q]
        qinn = self.mcmcpreverse[yy[q]]
        step_mat_big[pi,qi] = step_mat[p,q]
        
    self.step_mat = step_mat
    self.step_mat_big = step_mat_big
    return None


  def printstepmatrix(self,outfname=None):
    """
    print out the step matrix to outfname if not None, otherwise to [chainfname].step
    """
    if outfname is None:
      outfname = self.fname + '.step'

    nparam = len(self.mcmcp)
    ff = open(outfname,'w')
    for i in range(nparam):
      for j in range(nparam):
        ff.write('%15.5e' % self.step_mat_big[i,j])
      ff.write('\n')
    ff.close()

  def fillsig2hod(self,whichbox=1,cenopt=2,satopt=2,bosswdir='/home/howdiedoo/boss/',whichPB=0,COMVopt=0,cenvfrac=0.3,ihvscale=1.0,cenvfromchain=1,ihvfromchain=1):
    """
    Default behavior is to take values for cenvfrac and ihv from the chain, but add cenvfrac from the parameter in quadruture.
    """
    sig2hod = np.zeros(self.nelts) 
    mycenvfrac = cenvfrac
    myihv = ihvscale

    for cnt, elt in zip(range(self.nelts), self.chain):
      hh = hod.hodfromchain(elt,whichbox=whichbox,cenopt=cenopt,satopt=satopt,bosswdir=bosswdir,whichPB=whichPB,COMVopt=COMVopt) 
      if cenvfromchain == 1:
        mycenvfrac = (elt['cenvfrac']**2 + cenvfrac**2)
        if mycenvfrac > 0.:
          mycenvfrac = mycenvfrac**0.5
      if ihvfromchain == 1:
        myihv = elt['ihvscale']
      sig2hod[cnt] = hh.sig2hod(cenvfrac=mycenvfrac,ihvscale=myihv)
      #if cnt%100 == 0:
      #  print 'done',cnt,'of ',self.nelts
    self.sig2hod = sig2hod
    self.cenvfrac = cenvfrac
    self.ihvscale = ihvscale
    self.cenvfromchain = cenvfromchain
    self.ihvfromchain = ihvfromchain

  def writesig2hod(self,chainfname):
    outfname = chainfname + '.sig2hod'
    try:
      tmp = self.sig2hod[0]
    except:
      return 1
    if len(self.sig2hod) != len(self.chain['weight'][:]):
      return 1

    ofp = open(outfname,'w')
    ofp.write('# %d %d\n' % (self.cenvfromchain, self.ihvfromchain))
    ofp.write('# cenvfrac = %e\n' % (self.cenvfrac))
    if self.ihvfromchain != 1:
      ofp.write('# ihvscale = %e\n' % (self.ihvscale))
    else:
      if (np.fabs(self.chain['ihvscale'][:] - self.chain['ihvscale'][:].mean()) < 2.0e-5).all():
        ofp.write('# ihvscale = %e\n' % (self.chain['ihvscale'][:].mean()))
      else:
        ofp.write('# ihvscale = -1.0\n')
    for i in range(len(self.sig2hod)):
      ofp.write('%d %e\n' % (self.chain['weight'][i], self.sig2hod[i]))
    ofp.close()
    return 0

  def randomsubsample(self):
    """
    Returns a random element of the chain, where the probability to return an elt is weighted by weight.
    """
    wgtsum = self.chain['weight'][:].sum()
    rr = np.random.randint(0,wgtsum)
    ielt = np.where(self.chain['weight'][:].cumsum() >= rr)[0].min()
    return self.chain[ielt]


  def hodsubsample(self,marray,flist=[0.682689492137,0.954499736104],nsubsample=1000,
        whichbox=1,cenopt=2,satopt=2,bosswdir='/home/howdiedoo/boss/',whichPB=0,COMVopt=0):

    """
    Subsamples the chain and computes a bands containing flist fractions.
    nsubsample = 1000 by default.
    marray are mass points where you want to evaluate the hod.
    """
#def wgtresult2(wgt,v,flist=[0.682689492137,0.954499736104]):
    sNcen = np.zeros([len(marray),nsubsample])
    sNsat = np.zeros([len(marray),nsubsample])
    sNtot = np.zeros([len(marray),nsubsample])
    sMmed = np.zeros([3,nsubsample])
    sMavg = np.zeros([3,nsubsample])

    for nn in range(nsubsample):
      hh = hod.hodfromchain(self.randomsubsample(),whichbox=whichbox,cenopt=cenopt,satopt=satopt,\
                            bosswdir=bosswdir,whichPB=whichPB,COMVopt=COMVopt)
      sNcen[:,nn] = hh.Ncen(marray)
      sNsat[:,nn] = hh.Ncen(marray) * hh.Nsat(marray) 
      sNtot[:,nn] = sNcen[:,nn] + sNsat[:,nn] 
      ## compute interesting stats on halo masses.
      ctmp = (hh.mf.Nofm*hh.Ncen(hh.mf.m)).cumsum()
      sMmed[0,nn] = hh.mf.m[np.where(ctmp >= 0.5*ctmp[-1])[0].min()]
      sMavg[0,nn] = (hh.mf.m*hh.mf.Nofm*hh.Ncen(hh.mf.m)).sum()/ctmp[-1]
      ctmp = (hh.mf.Nofm*hh.Ncen(hh.mf.m)*(1.+hh.Nsat(hh.mf.m))).cumsum()
      sMmed[2,nn] = hh.mf.m[np.where(ctmp >= 0.5*ctmp[-1])[0].min()]
      sMavg[2,nn] = (hh.mf.m*hh.mf.Nofm*hh.Ncen(hh.mf.m)*(1.+hh.Nsat(hh.mf.m))).sum()/ctmp[-1]
      ctmp = (hh.mf.Nofm*hh.Ncen(hh.mf.m)*(hh.Nsat(hh.mf.m))).cumsum()
      sMmed[1,nn] = hh.mf.m[np.where(ctmp >= 0.5*ctmp[-1])[0].min()]
      sMavg[1,nn] = (hh.mf.m*hh.mf.Nofm*hh.Ncen(hh.mf.m)*(hh.Nsat(hh.mf.m))).sum()/ctmp[-1]


    rr = np.zeros([len(flist)*2+1, len(marray)])
    wgttmp = np.zeros(nsubsample)+1.
    bNcen = np.zeros([5,len(marray)])
    bNsat = np.zeros([5,len(marray)])
    bNtot = np.zeros([5,len(marray)])
 
    Mmed = np.zeros([5,3])
    Mavg = np.zeros([5,3])      
    for ii in range(3):    
      Mmed[:,ii] = wgtresult2raw(wgttmp,sMmed[ii,:])
      Mavg[:,ii] = wgtresult2raw(wgttmp,sMavg[ii,:])

    for mm in range(len(marray)):
#       print bNcen[:,mm].shape
#       print wgtresult2raw(wgttmp,sNcen[mm,:])
#       print wgtresult2raw(wgttmp,sNcen[mm,:]).shape
       bNcen[:,mm] = wgtresult2raw(wgttmp,sNcen[mm,:])
       bNsat[:,mm] = wgtresult2raw(wgttmp,sNsat[mm,:])
       bNtot[:,mm] = wgtresult2raw(wgttmp,sNtot[mm,:])
    self.bNcen = bNcen
    self.bNsat = bNsat
    self.bNtot = bNtot
    self.msub = marray
    self.Mmed = Mmed
    self.Mavg = Mavg
    self.nsubsample = nsubsample

  def plothodsubsample(self,ax=None,popt=0,rebinsize=10,fillopt=0,color='k',renormalize=False,MFwgt=False,flist=[0.954499736104],nsubsample=1000,
        whichbox=1,cenopt=2,satopt=2,bosswdir='/home/howdiedoo/boss/',whichPB=0,COMVopt=0,alpha=0.5,fsatrescale=1.):
    """
    Computes if necessary the subsample.
    popt = 0: Ncen
    popt = 1: Nsat
    popt = 2: Ntot
    popt = 3: Nsat/fsatrescale; only makes sense on a MFwgt plot.
    """
    try:
      x = self.bNcen[0,0]
      x = self.bNcen[0,0]
      x = self.bNcen[0,0]
      wgt = np.zeros(len(self.msub))+1.
      if MFwgt == True:
        hh = hod.hodfromchain(self.chain[0],whichbox=whichbox,cenopt=cenopt,satopt=satopt,bosswdir=bosswdir,whichPB=whichPB,COMVopt=COMVopt) 
        wgt = hh.mf.Nofm
        if rebinsize > 1:
          Nrebin = BRutils.rebin(hh.mf.Nofm,rebinsize=rebinsize)
          wgt = Nrebin
          
    except:
      ## need to run subsample.
      hh = hod.hodfromchain(self.chain[0],whichbox=whichbox,cenopt=cenopt,satopt=satopt,bosswdir=bosswdir,whichPB=whichPB,COMVopt=COMVopt) 
      Mlist = hh.mf.m
      wgt = np.zeros(len(Mlist))+1.
      if MFwgt == True:
        wgt = hh.mf.Nofm
      if rebinsize > 1:
        Mlist = np.exp(BRutils.rebin(np.log(hh.mf.m),rebinsize=rebinsize))
        wgt = np.zeros(len(Mlist))+1.
        tmpo = np.log10(Mlist[1:]/Mlist[:-1]).mean()
        assert (np.fabs(np.log10(Mlist[1:]) - np.log10(Mlist[:-1]) - tmpo) < 0.001*tmpo).all()
        if MFwgt == True:
          Nrebin = BRutils.rebin(hh.mf.Nofm,rebinsize=rebinsize)
          wgt = Nrebin
      self.hodsubsample(Mlist,flist=flist,nsubsample=nsubsample,whichbox=whichbox,cenopt=cenopt,satopt=satopt,bosswdir=bosswdir,whichPB=whichPB,COMVopt=COMVopt)
    if ax is None:
      ff = plt.figure(figsize=[6,6])
      ax = ff.add_subplot(1,1,1)
    else:
      ff = None
    if popt == 0:
      myy1 = wgt*self.bNcen[1]
      myy2 = wgt*self.bNcen[2]
    if popt == 1:
      myy1 = wgt*self.bNsat[1]
      myy2 = wgt*self.bNsat[2]
    if popt == 2:
      myy1 = wgt*self.bNtot[1]
      myy2 = wgt*self.bNtot[2]
    if popt == 3:
      myy1 = wgt*self.bNsat[1]/fsatrescale
      myy2 = wgt*self.bNsat[2]/fsatrescale

    if renormalize == True:
      myy1 = myy1/float(myy1.sum())
      myy2 = myy2/float(myy2.sum())

    if fillopt == 1:
      ax.fill_between(self.msub,myy1,myy2,facecolor=color,alpha=alpha)
    else:
      ax.plot(self.msub,myy1,color)
      ax.plot(self.msub,myy2,color)
## how do we fill between two curves instead.
    ax.set_xscale('log')
    return ff, ax
 


def combinesteps(m1,m2,m1rescale,m2rescale):
  """
  The purpose of this function is to take some partial step matrices and combine them to get a 
  first hack at varying all the parameters.
  This is a simple sum m1*m1rescale + m2*m2rescale.
  You have to account for changes in parameter numbers yourself in the rescale inputs.
  """
  if m1.shape != m2.shape:
    print 'matrices not aligned.',m1.shape,m2.shape
    return None 

  return m1*m1rescale + m2*m2rescale


def tablefmt1(ccx,myfs8,startstring=' 1 & Y & N '):
  mystr = '' + startstring
  for ndigits, nn in zip([3,2,2,2,2,2,4,3,3],['M_min','sigma_logM','M_cut', 'M1', 'alpha', 'nbar', 'fsat', 'hvscale', 'ihvscale', 'cenvfrac']):
    m,s = wgtresultint(ccx.chain['weight'][:], ccx.chain[nn][:])
    if re.match(nn,'hvscale'):
      m = m*myfs8 #= 0.48 ## need to compute/confirm this!!!
      s = s*myfs8
    if re.match(nn,'nbar'):
      m = 1.e4*m
      s = 1.e4*s
    #myfmt = '%.' + str(ndigits) + 'f \pm %.' + str(ndigits)+'f'
    #mystr = mystr + '& $'+str(myfmt % (m,s))+'$'
    myfmt = '%.' + str(ndigits)+'f'
    mystr = mystr + '& $ '+str(myfmt % m) + r' \pm ' + str(myfmt % s)+' $'
  mystr = mystr + r'\\'
  return mystr
    
def writetablefmt2(clist,clabel,fs8list,fname,boldfixedopt=0,boldcollist=[]):
  ofp = open(fname,'w')
  ofp.write(' ') ## no column name for first column of params.
  for cnt, xx in zip(range(len(clabel)), clabel):
    if not (cnt in boldcollist):
      ofp.write(r' & %s' % str(xx))
    else:
      ofp.write(r' & {\bf %s}' % str(xx))
  ofp.write('\\\\\n')
  for ndigits, nn, nnshow in zip([3,2,2,2,2,2,4,3,2,2,1,1,1],['M_min','sigma_logM','M_cut', 'M1', 'alpha', 'nbar', 'fsat', 'hvscale', 'ihvscale', 'cenvfrac','chi2_wp','chi2_multi','chi2_tot'],\
['$\log_{10} M_{\\rm min}$','$\sigma_{\log_{10} M}$','$\log_{10} M_{\\rm cut}$', '$\log_{10} M_1$', '$\\alpha$', '$\\bar{n}$', '$f_{\\rm sat}$', '$f\sigma_8$', '$\gamma_{\\rm IHV}$', '$\gamma_{\\rm cenv}$','$\\chi^2_{w_p}$ (18)', '$\\chi^2_{\\hat{\\xi}_{0,2}}$ (18)','$\\chi^2_{w_p+\\hat{\\xi}_{0,2}}$ (27)']):
    mystr = nnshow
    for cnt, ccx, myfs8, ccxl in zip(range(len(clabel)), clist,fs8list,clabel):
      m,s = wgtresultint(ccx.chain['weight'][:], ccx.chain[nn][:])
      if re.match(nn,'hvscale'):
        m = m*myfs8 #= 0.48 ## need to compute/confirm this!!!
        s = s*myfs8
      if re.match(nn,'nbar'):
        m = 1.e4*m
        s = 1.e4*s
    
      myfmt = '%.' + str(ndigits)+'f'
      if re.search('chi2',nn):
        xxx = np.where(ccx.chain['chi2_tot'][:] == ccx.chain['chi2_tot'][:].min())[0]
        m = ccx.chain[nn][xxx[0]]
        if cnt in boldcollist:
          mystr = mystr + r' & ${\bf '+str(myfmt % m) +' }$'
        else:
          mystr = mystr + ' & $ '+str(myfmt % m) +' $'
        #print 'hello beth!'
      elif nn not in ccx.mcmcp.keys() or ccx.mcmcfixed[ccx.mcmcp[nn]] == 0:
        if cnt in boldcollist:
          mystr = mystr + r' & ${\bf '+str(myfmt % m) + r' \pm ' + str(myfmt % s)+' }$'
        else:
          mystr = mystr + ' & $ '+str(myfmt % m) + r' \pm ' + str(myfmt % s)+' $'

      else:
        if boldfixedopt == 0 and not (cnt in boldcollist):
          mystr = mystr + ' & $ '+str(myfmt % m) +' $'
        else:
          mystr = mystr + r' & ${\bf '+str(myfmt % m) + r' }$'

      
    mystr = mystr + r'\\'
    ofp.write("%s\n" % (mystr))
  ofp.close()
  #return mystr

