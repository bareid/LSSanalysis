import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ximisc

class xiell:
  def __init__(self,xiellfname=None,icovfname=None,ellmax=-1,sxilist=[],smincut=-1.,smaxcut=1.e12,skiprows=0,srescale=1.):
    """
    If you only want to read in some of the ell values, use ellmax.
    Otherwise reads in all the multipoles in the file.
    Use smincut/smaxcut to remove outlying bins you dont' want to keep.
    srescale is to convert theory files in Mpc to Mpc/h.  Rescale value is h (0.7)"
    """
#    try:
    if(0==0):
      if(xiellfname is not None):
        self.fname = [xiellfname]
        ## accept 2 formats.
  ## later extend format to include xi0,2,4?
        aa = np.loadtxt(xiellfname,skiprows=skiprows)
        if(np.fabs(srescale-1.) > 1.0e-5):
          aa[:,0] = aa[:,0]*srescale
        s = aa[:,0]
        ikeep = np.where((s >= smincut) & (s <= smaxcut))[0]
        aa = aa[ikeep,:]
        s = aa[:,0]
        ns = len(s)
        nxi = len(aa[0,:])-1
        assert nxi > 0
        if nxi == 1:  # long format, s, xi0 followed by s, xi2.
          splitz = np.where(s[:-1] > s[1:])[0]
          nell = len(splitz) + 1
          if(len(splitz) == 0):
            nlist = [len(s)]
          else:
            nlist = []
            ilow = 0
            for sp in splitz:
              ihigh = sp+1
              nlist.append(ihigh-ilow)
              ilow=ihigh
            nlist.append(len(s)-ilow)
          xi=np.zeros([nell,max(nlist)])
          svec=np.zeros([nell,max(nlist)])
          ilow = 0
          elli = 0
          for n in nlist:
            ihigh = ilow + n
            svec[elli,:n] = aa[ilow:ihigh,0]
            xi[elli,:n] = aa[ilow:ihigh,1]
            if(elli == 0):
              self.s0 = svec[elli,:n]
              self.xi0 = xi[elli,:n]
            if(elli == 1):
              self.s2 = svec[elli,:n]
              self.xi2 = xi[elli,:n]
            ilow = ihigh
            elli += 1
  
  ## use xilong to compute chi2
        self.slong = aa[:,0]
        self.xilong = aa[:,1]
        self.ndata = len(aa[:,1])
  
  
        if(nxi > 1):
          nell = nxi
          if(ellmax != -1 and (ellmax == 0 or ellmax == 2 or ellmax == 4)):
            nell = ellmax/2+1
          
          xi = np.zeros([nell,ns])
          svec = np.zeros([nell,ns])
          xilong = []
          self.slong = aa[:,0]
       
          for elli in range(nell):
            svec[elli,:] = aa[:,0]
            xi[elli,:] = aa[:,elli+1]
            xilong.append(aa[:,elli+1])
            if(elli == 0):
              self.s0 = svec[elli,:]
              self.xi0 = xi[elli,:]
            if(elli == 1):
              self.s2 = svec[elli,:]
              self.xi2 = xi[elli,:]
  
          xilong = np.array(xilong).flatten()
          self.xilong = xilong
          self.ndata = len(self.xilong)
  
        self.svec = svec
        self.xi = xi
        self.nell = nell

### create a xiell object from numpy arrays.
      elif sxilist != []:
        if(np.fabs(srescale-1.) > 1.0e-5):
          sxilist[0] = sxilist[0]*srescale
        self.svec = sxilist[0]
        self.xi = sxilist[1]
        self.nell = len(self.svec[:,0])
        self.ndata = len(self.xi.flatten())
        self.xilong = self.xi.flatten()
        self.fname = None
        self.s0 = self.svec[0,:]
        self.xi0 = self.xi[0,:]
        self.s2 = self.svec[2,:]
        self.xi2 = self.xi[2,:]

      else:
        assert 0==1


#    except:
    else:
      print 'bad xiell file',xiellfname
      self = None

    if(self is not None):
      ## let's finde out if this is a data of theory curve.
      try:
        fbase = xiellfname.split('.'+xiellfname.split('.')[-1])[0]
        ifpDR = open(fbase+'.DRopt2','r')
        line = ifpDR.readline()
        line = ifpDR.readline()
        DDwgt = float(line.split(':')[1].split(',')[0])
        RRwgt = float(line.split(':')[1].split(',')[1])
        DRfac = float(DDwgt)/float(RRwgt)
        ifpDR.close()
        ## set weight to sum of data weights.
        self.weight = DDwgt
## didn't work, not sure why.
      #except IOerror as e:
      except:
#        print 'could not open',fbase+'.DRopt2'
#        print 'so setting weight to 1.0'
        self.weight = 1.0

    if(self is not None):
      if(icovfname is None):
        self.DorT = 1
      else:
        self.DorT = 0
        icov = np.matrix(np.loadtxt(icovfname))
        self.icov = icov
        assert len(icov[:,0]) == self.ndata
## compute diag errors.
        cov = icov.I
        diagerr = np.array([np.sqrt(cov[i,i]) for i in range(self.ndata)])
        self.diagerr = diagerr
        self.diagerrxi0 = diagerr[0:len(self.xi0)]
        self.diagerrxi2 = diagerr[len(self.xi0):len(self.xi0)+len(self.xi2)]


  def __str__(self):
    mystr = 'xiell (nell = %d) data vector from %s\n' % (self.nell, self.fname)
    for i in range(self.ndata):
      mystr = mystr +  '%e %e\n' % ((self.svec.ravel())[i],self.xilong[i])
    print mystr
    mystr = 'hi'
    return mystr

  

  def checkbinning(self,other):
    """
    Returns 0 if binning is the same, 1 otherwise
    """
    if(self.ndata != other.ndata):
      return 1
    if (np.fabs(self.svec.flatten()[0:self.ndata] - other.svec.flatten()[0:self.ndata]) > 2.0e-4).any():
      return 1
    return 0

  def checkbinnings02(self,other):
    """
    Returns 0 if binning is the same in xi0 and xi2, 1 otherwise.
    """
    if(len(self.s0) != len(other.s0)):
      return 1
    if(len(self.s2) != len(other.s2)):
      return 1
    if(np.fabs(self.s0[:] - other.s0[:]) > 2.0e-4).any():
      return 1 
    if(np.fabs(self.s2[:] - other.s2[:]) > 2.0e-4).any():
      return 1 
    return 0

  def __add__(self,other):
    if(self.checkbinning(other) == 0):
      xiellnew = copy.deepcopy(self)
      xiellnew.xi = (self.xi*self.weight + other.xi*other.weight)/(self.weight+other.weight)
      xiellnew.xilong = (self.xilong*self.weight + other.xilong*other.weight)/(self.weight+other.weight)
      xiellnew.weight = (self.weight+other.weight)

    else:
      print 'add fail!  binning misaligned'
      xiellnew = None
    return xiellnew

  def __sub__(self,other):
    if(self.checkbinning(other) == 0):
      xiellnew = copy.deepcopy(self)
      xiellnew.xi = other.xi-self.xi
      xiellnew.xilong = other.xilong-self.xilong
    else:
      print 'subtract fail!  binning misaligned'
      xiellnew = None
    return xiellnew


  def chi2(self,other):
    """
    compute chi2 of self with other.  self is the one with the inverse cov.
    """
    diff = np.matrix(self.xilong-other.xilong)
    if self.DorT == 0 and other.DorT == 1:
      return ((diff*self.icov) * (diff.T))[0,0]
    elif self.DorT == 1 and other.DorT == 0:
      return ((diff*other.icov) * (diff.T))[0,0]
    else:
      return 0. ## can't compute chi2 without an inverse covariance matrix.

  def printxiellshort(self,outfname):
          #svec = np.zeros([nell,ns])
    ## first check whether this print makes sense -- all the svec values are the same.
    for elli in range(self.nell): 
      for ellj in range(elli+1,self.nell): 
        if (np.fabs(self.svec[elli,:] - self.svec[ellj,:]) > 2.0e-5).any():
          return None 
    ofp = open(outfname,'w')
    for i in range(len(self.svec[0,:])):
      ofp.write('%e' % (self.svec[0,i]))
      for ell in range(self.nell):
        ofp.write(' %e' % (self.xi[ell,i]))
      ofp.write('\n')
    ofp.close()

  def printxielllong(self,outfname):
    for i in range(self.ndata):
      ofp.write('%e %e\n' % (self.slong[i], self.xilong[i])) 

  def addcurve(self,ax,ell=0,color='k',spow=1,fmt=None,lbl=None):
    """
    adds a curve on a xiell plot. specify which ell you want with ell.
    """
    stmp = self.svec[ell/2,:]
    if(self.DorT == 0 and (ell == 0 or ell == 2)): ## data, so there are errors.
      if(ell == 0):
        err = self.diagerrxi0
      if(ell == 2):
        err = self.diagerrxi2
      ax.errorbar(stmp,stmp**spow*self.xi[ell/2,:],yerr=err*stmp**spow,color=color,ecolor=color,fmt=fmt,label=lbl)
    else:
      ax.plot(stmp,stmp**spow*self.xi[ell/2,:],label=lbl,color=color)

  def makeplot(self,ell=0,ax=None,color='k',span=None,logxopt=1,logyopt=0,spow=1,fmt=None,lbl=None):
    """
    Makes a standard xiell 1d plot.
    log[x-y]opt specify if you want log or linear binning on the axis.
    spow means multiply y coordinate by s**spow.
    """

    if ax is None:
      ## this runs into some self-referential problems when copying class objects, so taking this out.
      #self.fig = plt.figure(figsize=[6,6])
      ff = plt.figure(figsize=[6,6])
      ax=ff.add_subplot(1,1,1)
    else:
      ff = None

    self.addcurve(ax,ell,color=color,spow=spow,lbl=lbl,fmt=fmt)

    if span is None:
      span = [self.svec[ell/2,:].min()*0.9,self.svec[ell/2,:].max()*1.1,(self.svec[ell/2,:]**spow*self.xi[ell/2,:]).min()*0.9, (self.svec[ell/2,:]**spow*self.xi[ell/2,:]).max()*1.1]
    if(logxopt == 1):
      ax.set_xscale('log')
    if(logyopt == 1):
      ax.set_yscale('log')


    ax.axis(span)
    ax.set_xlabel(r'$s \, [h^{-1} {\rm Mpc}]$',fontsize=20)
    if(np.fabs(spow) > 0.01):
      ax.set_ylabel(r'$s^{%.1f} \xi_{%d}(s)$' % (spow, ell),fontsize=20)
    else:
      ax.set_ylabel(r'$\xi_{%d}(s)$' % (ell),fontsize=20)

    return ff, ax

  def testo(self,other):
    """debugging copy problems."""
    print 'i want to copy this guy:'
    print self
    a1=copy.deepcopy(self) 
    print 'i want to copy this guy:'
    print other
    a2=copy.deepcopy(self)

    print a1.xi[2,40]
    print a2.xi[2,40]

  def makediffplot(self,other,ell=0,ax=None,color='k',span=None,logxopt=1,logyopt=0,spow=1,fmt=None,lbl=None,ecolor='k'):
    """
    Makes a standard xiell 1d plot, but plot difference between data (self) and theory (other).
    log[x-y]opt specify if you want log or linear binning on the axis.
    spow means multiply y coordinate by s**spow.
    """

    if ax is None:
      ## this runs into some self-referential problems when copying class objects, so taking this out.
      #self.fig = plt.figure(figsize=[6,6])
      ff = plt.figure(figsize=[6,6])
      ax=ff.add_subplot(1,1,1)
    else:
      ff = None

    ## want to plot the errors about 0, use a tmp
    tmp1 = copy.deepcopy(self)
    tmp1.xi[:,:] = 0.
    tmp1.xilong[:] = 0.

    tmp1.addcurve(ax,ell,color=ecolor,spow=spow,lbl=None,fmt=fmt)

    tmp2 = self - other
    tmp2.DorT = 1
    tmp2.addcurve(ax,ell,color=color,spow=spow,lbl=lbl,fmt=fmt)

    ymin = min((self.svec[ell/2,:]**spow*tmp1.xi[ell/2,:]).min()*0.9, \
              (self.svec[ell/2,:]**spow*tmp2.xi[ell/2,:]).min()*0.9)
    ymax = max((self.svec[ell/2,:]**spow*tmp1.xi[ell/2,:]).max()*1.1, \
              (self.svec[ell/2,:]**spow*tmp2.xi[ell/2,:]).max()*1.1)


    if span is None:
      span = [self.svec[ell/2,:].min()*0.9,self.svec[ell/2,:].max()*1.1,ymin,ymax]
    if(logxopt == 1):
      ax.set_xscale('log')
    if(logyopt == 1):
      ax.set_yscale('log')


    ax.axis(span)
    ax.set_xlabel(r'$s \, [h^{-1} {\rm Mpc}]$',fontsize=20)
    if(np.fabs(spow) > 0.01):
      ax.set_ylabel(r'$s^{%.1f} \Delta \xi_{%d}(s)$' % (spow,ell),fontsize=20)
    else:
      ax.set_ylabel(r'$\Delta \xi_{%d}(s)$' % (ell),fontsize=20)

    del tmp1
    del tmp2
    return ff, ax

  #### proper indentation -- cannot be part of the class.
  ## copying from boss/zdistvXlogbinsompcleverLSnew/xitruegenericallscalesRRdown.py for this function.
  ## and also boss/zdistvXlogbinsompcleverLSsmallscale/rebinxismugenericcutsmallscale.py

def xiellfromDR(fbase,nell=3,binfile=None,rperpcut=-1.,dfacs=1,dfacmu=1,icovfname=None,smincut=-1.,smaxcut=1.e12):
  """
  This function should supercede older codes with names like rebinxismugeneric*py
  input the base for the DR files, as well as optional bin files and small scale 
  cutoff for masking (rperpcut).  If you don't want to use a binfile to specify
  the downsampling, use dfacs and dfacmu inputs.
  smincut/smaxcut will only be used if binfile is None.
  """
  fDD = fbase+'.DRopt1'
  fDR = fbase+'.DRopt2'
  fRR = fbase+'.DRopt3'
  rg, mug, DDg = np.loadtxt(fDD,skiprows=3,unpack=True)
  rg, mug, DRg = np.loadtxt(fDR,skiprows=3,unpack=True)
  rg, mug, RRg = np.loadtxt(fRR,skiprows=3,unpack=True)
  if binfile is None:
    ikeep = np.where((rg >= smincut) & (rg <= smaxcut))[0]
    rg = rg[ikeep]
    mug = mug[ikeep]
    DDg = DDg[ikeep]
    DRg = DRg[ikeep]
    RRg = RRg[ikeep]
  ifpRR = open(fRR,'r')
  line = ifpRR.readline()
  line = ifpRR.readline()
  ## rerun everything to get this to 12 digits, just in caes!
  DDwgt = float(line.split(':')[1].split(',')[0])
  RRwgt = float(line.split(':')[1].split(',')[1])
  RRfacdown = float(RRwgt)
  ifpRR.close()

  ifpDR = open(fDR,'r')
  line = ifpDR.readline()
  line = ifpDR.readline()
  ## rerun everything to get this to 12 digits, just in caes!
  DDwgt = float(line.split(':')[1].split(',')[0])
  RRwgt = float(line.split(':')[1].split(',')[1])
  DRfac = float(DDwgt)/float(RRwgt)
  #print DRfac
  ifpDR.close()

  fixRR = RRwgt/RRfacdown
  #print 'fixRR: ',fixRR

  nmu = len(np.where(rg == rg[0])[0])
  dmu  = 1./float(nmu)
  nr = len(np.where(mug == mug[0])[0])
  if nmu*nr != len(rg):
    return None

  ## is s (called r here) log or linear binning?
  s1d = rg.reshape([nr,nmu])[:,0]
  mu1d = rg.reshape([nr,nmu])[0,:]
  ## compute linear spacings.
  dmutmp = (mu1d[1:]-mu1d[:-1]).mean()
  ds = (s1d[1:]-s1d[:-1]).mean()
  if (np.fabs(mu1d[1:]-mu1d[:-1] - dmutmp) > 0.0001*dmutmp).any():
    ## not linear binning!
    ## this shouldn't be the case for a s,mu grid.
    assert 0==1
    assert np.fabs(dmutmp - dmu) < 1.0e-4

  if (np.fabs(s1d[1:]-s1d[:-1] - ds) > 0.0001*ds).any():
    logsopt = 0
    rglowedge = rg.copy()
    rglowedge = rglowedge-0.5*ds
  dlogs = (np.log(s1d[1:]/s1d[:-1])).mean()
  if (np.fabs(np.log(s1d[1:])-np.log(s1d[:-1]) - dlogs) < 0.0001*dlogs).all():
    logsopt = 1
    rglowedge = rglowedge*np.exp(-0.5*dlogs)

  assert logsopt == 0 or logsopt == 1

  ## cut at edge of bin.
  xx = np.where(rglowedge*(1-(mug+0.5*dmu)**2)**0.5 < rperpcut)[0]
  if(rperpcut < 0.):  assert len(xx) == 0
  if len(xx) > 0:
    print 'hi beth, cutting this many bins:',len(xx)
    DDg[xx] = 0.
    DRg[xx] = 0.
    RRg[xx] = 0.

  if binfile is not None:
    ## specifies how many rbins to join together for first bin, next bin, etc.
    rjoin, mujoin = np.loadtxt(binfile,unpack=True,dtype='int')
    nrcut = len(rjoin)
    ### use later.
    #bintag = binfile.split('.')[0]
  else:
#    dfacs = dfacsin
#    dfacmu = dfacmuin
    nrcut = len(s1d)/dfacs
    dmudown = dfacmu*dmu

  xi = np.zeros([nell,nrcut],dtype=float)
  svec = np.zeros([nell,nrcut],dtype=float)

  rcen = np.zeros(nrcut,dtype=float)
  ristart = 0
  for i in range(nrcut):
    if binfile is not None:
      dfacs = rjoin[i]
      dfacmu = mujoin[i]
      dmudown = dfacmu*dmu
    riend = ristart + dfacs - 1
    if logsopt == 0:
      rcen[i] = 0.5*(s1d[ristart] + s1d[riend])
    if logsopt == 1:
      rcen[i] = np.exp(0.5*np.log(s1d[ristart]*s1d[riend]))
    for ishort in range(dfacs):
      i1 = (ristart+ishort)*nmu
      i2 = (ristart+ishort+1)*nmu
      if(ishort == 0):
        mymu = ximisc.downsample1d(mug[i1:i2],dfacmu)
        mydd = ximisc.downsample1d(DDg[i1:i2],dfacmu)
        mydr = ximisc.downsample1d(DRg[i1:i2],dfacmu)
        myrr = ximisc.downsample1d(RRg[i1:i2],dfacmu)
      else:
        mydd = mydd + ximisc.downsample1d(DDg[i1:i2],dfacmu)
        mydr = mydr + ximisc.downsample1d(DRg[i1:i2],dfacmu)
        myrr = myrr + ximisc.downsample1d(RRg[i1:i2],dfacmu)

    yy = np.where(myrr < 0.01)[0]
    for ell in [0,2,4]:
      ###### WRONG!! #####
      #xi[i,ell/2] = ((mydd-mydr*DRfac)/myrr/DRfac**2*ximisc.legendre(ell,mymu)).sum()*dmudown*(2.*ell+1.)
      ### correct, but can't do fixRR
#        xitmp = (mydd-mydr*DRfac)/myrr/DRfac**2+1.
      xitmp = (mydd-mydr*DRfac)/myrr/DRfac**2/fixRR**2+1.
      xitmp[yy] = 0.  ## fix 0'd out regions
      xi[ell/2,i] = (xitmp*ximisc.legendre(ell,mymu)).sum()*dmudown*(2.*ell+1.)
      svec[ell/2,i] = rcen[i]
    ristart = ristart + dfacs

  return xiell(sxilist=[svec, xi],icovfname=icovfname)


if __name__ == '__main__':
  print 'hi'

  binfile="/Users/bareid/work/montserratdata/boss/zdistvXlogbinsompcleverLSnew/bin1match.txt"
  ival=1
  d2="/Users/bareid/work/montserratdata/boss/zdistvXlogbinsompcleverLSnew/outputQPMmocks/"
  fbase=d2+"a0.6452_%04d.dr10_ngc_rmax160deltalog10r" % ival
  tt=xiellfromDR(fbase=fbase,binfile=binfile)
  tt.printxiellshort("testo.dat")

  print 'exists?',tt.checkbinnings02(tt)
 
  ## this passed a test against the old version of the code::
#bash-3.2$ diff a0.6452_0001.dr10_ngc_rmax160deltalog10rrebin-bin1match.xiellcut ~/gitprojects/LSSanalysis/testo.dat 
#4c4
#< 2.851019e+00 4.957053e+00 4.155953e+00 4.444305e+00
#---
#> 2.851018e+00 4.957053e+00 4.155953e+00 4.444305e+00
#9c9
#< 2.317395e+01 1.964722e-01 -1.200379e-01 1.204180e-02
#---
#> 2.317394e+01 1.964722e-01 -1.200379e-01 1.204180e-02
  #



