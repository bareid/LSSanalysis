import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ximisc
import cov

class xiell:
  def __init__(self,xiellfname=None,icovfname=None,covfname=None,ellmax=-1,sxilist=[],smincut=-1.,smaxcut=1.e12,skiprows=0,srescale=1.,unbiasicovfac=1.,unbiascovfac=1.):
    """
    If you only want to read in some of the ell values, use ellmax.
    Otherwise reads in all the multipoles in the file.
    Use smincut/smaxcut to remove outlying bins you don't want to keep.
    srescale is to convert theory files in Mpc to Mpc/h.  Rescale value is h (0.7)"
    unbiasicovfac (if not 1., float(ntot - p - 2)/float(ntot-1) from 0608064).
    remember if you do that then cov and icov won't be inverses of each other!!
    You can either do unbiasicovfac in your icov file (so no correction here)
    but then your diagerrs might be wrong (but your chi2 will be correct).
    Perhaps apply the inverse as a factor to cov?
    """
#    print 'beth is it a problem that this module imports cov and has a local variable named cov?  should i change it?'
#    try:
    if(0==0):
      if((xiellfname is not None) and (sxilist == [])):
        self.fname = [xiellfname]
        ## accept 2 formats.
  ## later extend format to include xi0,2,4?
        aa = np.loadtxt(xiellfname,skiprows=skiprows)
        if(np.fabs(srescale-1.) > 1.0e-5):
          aa[:,0] = aa[:,0]*srescale
        s = aa[:,0]
        stmp = (aa[:,0]).copy()
        stmp1 = (aa[:,0]).copy()
        ikeep = np.where((s >= smincut) & (s <= smaxcut))[0]
        nstart = len(s)
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
          stmp = np.zeros(nell*len(stmp1))  ## needed for cov stuff later.
          nstart = len(stmp) ## account for s values for each multipole
          for iii in range(nell):
            stmp[iii*len(aa[:,0]):(iii+1)*len(stmp1)] = stmp1[:]
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
          self.slong = svec.flatten()
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
        self.slong = self.svec.flatten()
        self.fname = None
        self.s0 = self.svec[0,:]
        self.xi0 = self.xi[0,:]
        if(self.nell > 1):
          self.s2 = self.svec[1,:]
          self.xi2 = self.xi[1,:]

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
      if(icovfname is None and covfname is None):
        self.DorT = 1
      else:
        self.DorT = 0
        if(icovfname is not None):
          icov = np.matrix(np.loadtxt(icovfname))
          self.icov = icov
          cov = icov.I
          self.cov = cov
        else:
          assert (covfname is not None)
          ### deal with antonio cov file format.
          covtmp = np.loadtxt(covfname)
          xshape, yshape = covtmp.shape
          if(yshape == 3):
            if(xshape == 3):
              print 'from current code, this shape is ambiguous.  Fix the code before proceeding!'
              self = None
              sys.exit(1)
            cov = np.loadtxt(covfname,usecols=[2],unpack=True)
            if sxilist == []: 
              assert (len(cov) == (self.ndata)**2) | (len(cov) == nstart**2)
            else:
              assert (len(cov) == (self.ndata)**2)
            cov = np.matrix(cov.reshape(int((len(cov)**0.5)), int((len(cov)**0.5))))
          else: # just read in regular matrix.
            cov = np.matrix(np.loadtxt(covfname))
            if sxilist == []: 
              assert (len(cov) == (self.ndata)) | (len(cov) == nstart)
            else:
              assert (len(cov) == (self.ndata))
          self.cov = cov
          self.icov = cov.I
        if(len(cov[:,0]) > self.ndata):  ## need to select desired rows/cols from the covariance matrix.
          cov = np.array(cov)
#tmp!
          print self.ndata, len(cov[:,0])
          x, y = np.meshgrid(stmp,stmp)
          ikeep = np.where((x >= smincut) & (x <= smaxcut) & (y >= smincut) & (y <= smaxcut))
          assert (x[ikeep] >= smincut).all()
          assert (y[ikeep] >= smincut).all()
          assert (x[ikeep] <= smaxcut).all()
          assert (y[ikeep] <= smaxcut).all()
          assert len(cov[ikeep]) == (self.ndata)**2
          cov = cov[ikeep].reshape(self.ndata,self.ndata)
          cov = np.matrix(cov)
          self.cov = cov
          self.icov = cov.I
          assert len(self.cov[:,0]) == self.ndata 
##  apply unbiasicovfac.  remember now icov != inverse of cov !!!
        if np.fabs(unbiasicovfac - 1.) > 2.0e-6 and np.fabs(unbiascovfac - 1.) > 2.0e-6:
          print 'only one unbiasfac can be set.'
          self = None
          sys.exit(1)
        self.icov = self.icov*unbiasicovfac
        self.unbiasicovfac = unbiasicovfac
        self.cov = self.cov*unbiascovfac
        self.unbiascovfac = unbiascovfac

## compute diag errors.
        diagerr = np.array([np.sqrt(cov[i,i]) for i in range(self.ndata)])
        self.diagerr = diagerr
        self.diagerrxi0 = diagerr[0:len(self.xi0)]
        if(self.nell > 1):
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
    elif(self.ndata == other.ndata):
      print 'strong warning!  binnings do not agree.  Subtracting anyway'
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
    can also input a long vector and will compute chi2 with htat if format correct.
    """
    if isinstance(other,xiell):
      diff = np.matrix(self.xilong-other.xilong)
      if self.DorT == 0 and other.DorT == 1:
        return ((diff*self.icov) * (diff.T))[0,0]
      elif self.DorT == 1 and other.DorT == 0:
        return ((diff*other.icov) * (diff.T))[0,0]
      elif self.DorT == 0 and other.DorT == 0:
        return ((diff*self.icov) * (diff.T))[0,0]
      else:
        return 0. ## can't compute chi2 without an inverse covariance matrix.
    else:
      try:
        diff = np.matrix(self.xilong - other)
      except:
        return 0.
      return ((diff*self.icov) * (diff.T))[0,0]

  def chi2info(self,other):
    """
    For non-diagonal cov, this helps show you which differences between model and theory contribute the most
    to chi2.
    """
    if isinstance(other,xiell):
      diff = np.matrix(self.xilong-other.xilong)
      if self.DorT == 0 and other.DorT == 1:
        myicov = self.icov
#        return ((diff*self.icov) * (diff.T))[0,0]
      elif self.DorT == 1 and other.DorT == 0:
        myicov = other.icov
#        return ((diff*other.icov) * (diff.T))[0,0]
      elif self.DorT == 0 and other.DorT == 0:
        myicov = self.icov
#        return ((diff*self.icov) * (diff.T))[0,0]
      else:
        return 0. ## can't compute chi2 without an inverse covariance matrix.
    else:
      try:
        diff = np.matrix(self.xilong - other)
        myicov = self.icov
      except:
        return 0.
    mysum = 0.
    for i in range(len(self.xilong)):
      print i,diff[0,i]*((diff*myicov)[0,i])
      mysum += diff[0,i]*((diff*myicov)[0,i])
    assert np.fabs(mysum - self.chi2(other)) < 0.001
    return mysum



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
    ofp = open(outfname,'w')
    for i in range(self.ndata):
      ofp.write('%e %e\n' % (self.slong[i], self.xilong[i])) 
    ofp.close()

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

  def makecosmoxi2dinput(self,outfbase,hfid=0.7):
    logorlinear = -1
    ## check for log or linear.
    ds0 = (self.s0[1:] - self.s0[:-1]).mean()

    if (np.fabs(self.s0[1:]-self.s0[:-1]-ds0) < 2.0e-4).all():
      if (self.nell > 1):
        ds2 = (self.s2[1:] - self.s2[:-1]).mean()
        if (np.fabs(self.s2[1:]-self.s2[:-1]-ds2) < 2.0e-4).all():
          logorlinear = 1
      else:
        logorlinear = 1

    ## check for log or linear xi0 and xi2
    dlogs0 = (np.log(self.s0[1:]/self.s0[:-1])).mean()
    if (np.fabs(np.log(self.s0[1:]/self.s0[:-1])-dlogs0) < 2.0e-4).all():
      if (self.nell > 1):
        dlogs2 = (np.log(self.s2[1:]/self.s2[:-1])).mean()
        if (np.fabs(np.log(self.s2[1:]/self.s2[:-1])-dlogs2) < 2.0e-4).all():
          logorlinear = 0
      else:
        logorlinear = 0

    if logorlinear == -1:
      print 'wacky binning, make your own bin file by hand.'
      return None

    ds0 = ds0/hfid
    if(self.nell > 1):
      ds2 = ds2/hfid

    ofp = open(outfbase+".bin",'w')
    ofpd = open(outfbase+".dat",'w')
    ofp.write("# hfid = %.6e\n" % (hfid))

    for i in range(len(self.s0)):
      if logorlinear == 1:
        ofp.write('%.7e %.7e %.7e\n' % (self.s0[i]/hfid, self.s0[i]/hfid-0.5*ds0,  self.s0[i]/hfid+0.5*ds0))
      else:
        ofp.write('%.7e %.7e %.7e\n' % (self.s0[i]/hfid, self.s0[i]/hfid*np.exp(-0.5*dlogs0),  self.s0[i]/hfid*np.exp(0.5*dlogs0)))
      ofpd.write("%15.7e %15.7e " % (self.s0[i]/hfid, self.xi0[i]))
      if(self.nell > 1):
        ofpd.write("%16.7e " % (self.xi2[i]))
      else:
        ofpd.write("%16.7e\n" % (0.0))

    ofp.close()
    ofpd.close()
    if(self.nell == 1):
      return 0
    else:
      assert (np.fabs(self.s0 - self.s2) < 2.0e-4).all()
    return 0

  #### proper indentation -- cannot be part of the class.
  ## copying from boss/zdistvXlogbinsompcleverLSnew/xitruegenericallscalesRRdown.py for this function.
  ## and also boss/zdistvXlogbinsompcleverLSsmallscale/rebinxismugenericcutsmallscale.py

def weightsum(xiinlist,wgtlist):
  """
  returns a new xiell instance whose xiell values are the weighted sum of 
  the values in the xiinlist. 
  """
  if len(xiinlist) < 1:
    return None
  try:
    ximean = copy.deepcopy(xiinlist[0].xi)
  except:
    return None

  ximean = ximean*wgtlist[0]
  wgtsum = wgtlist[0]
  for ii in range(1,len(xiinlist)):
## first sanity checks that binning is equal.
    if xiinlist[ii].nell != xiinlist[0].nell:
      return None
    if xiinlist[ii].ndata != xiinlist[0].ndata:
      return None
    
    for ll in range(xiinlist[0].nell):
      if (xiinlist[ii].svec[ll,:] != xiinlist[0].svec[ll,:]).any():
        return None
      if ll == 0:
        if (xiinlist[ii].s0 != xiinlist[0].s0).any():
          return None
      if ll == 1:
        if (xiinlist[ii].s2 != xiinlist[0].s2).any():
          return None

    ximean += xiinlist[ii].xi*wgtlist[ii]
    wgtsum += wgtlist[ii]


  ximean = ximean/float(wgtsum)
  return xiell(sxilist=[xiinlist[0].svec, ximean])

def weightsumvar(xiinlist,wgtlist):
  """
  returns diagonal variance, correcting for n/(n-1) factor.
  """
  ximean = weightsum(xiinlist,wgtlist)
  ximeanchk = np.zeros_like(ximean.xilong)
  xivar = np.zeros_like(ximean.xilong)
  for ii in range(len(xiinlist)):
    ximeanchk[:] += xiinlist[ii].xilong[:]
    xivar[:] += xiinlist[ii].xilong[:]*xiinlist[ii].xilong[:]

  ximeanchk = ximeanchk/float(len(xiinlist))
  xivar = xivar/float(len(xiinlist))
  assert (np.fabs(ximeanchk - ximean.xilong) < 2.0e-6).all()

  xivar = xivar - ximeanchk**2
  return xivar*float(len(xiinlist))/float(len(xiinlist)-1.)


## has this function been tested??

def xiellfromDR(fbase,nell=3,binfile=None,rperpcut=-1.,
                dfacs=1,dfacmu=1,icovfname=None,smincut=-1.,smaxcut=1.e12,\
                DRfacinfo=None,smallRRcut=-1.,periodicopt=0,printmaskopt=0,xiopt=0):
  """
  This function should supercede older codes with names like rebinxismugeneric*py
  input the base for the DR files, as well as optional bin files and small scale 
  cutoff for masking (rperpcut).  If you don't want to use a binfile to specify
  the downsampling, use dfacs and dfacmu inputs.
  smincut/smaxcut will only be used if binfile is None.
  DRfacinfo is used to pass [DRfac, fixRR] (dont get it from file headers in this case)
  That feature is necessary for computing covariance matrices from bootstrap, where that factor
  should be fixed by total N and/or S, it does not vary from bootstrap region to region.
  smallRRcut allows you to replace randoms(mu) with the average over mu where the total number of 
  randoms in the bin is smaller than smallRRcut, just to keep Poisson noise down.
  If periodicopt == 1, assume the input is file containing [rg, mug, Vg, Ng] as currently output
  by correlationfxnMASTERv4.
  xiopt introduced on May 13 2014 to input xi(s,mu) and simply rebin (for comparison with Hong).
  """
  if periodicopt == 1: # periodic box.
    rg, mug, Vg, Ng = np.loadtxt(fbase,unpack=True)
    Vzerochk = Vg.min()*0.1
    fDR = fbase ## just for the return xiell purposes.
  elif xiopt == 1:
    rg, mug, xig = np.loadtxt(fbase,unpack=True,usecols=[0,1,2])  ## Hong format.
    fDR = fbase
    if binfile is None:
      ikeep = np.where((rg >= smincut) & (rg <= smaxcut))[0]
      rg = rg[ikeep]
      mug = mug[ikeep]
      xig = xig[ikeep]
  else: # DD,RR counts.
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

    if DRfacinfo is None:
      DRfac, fixRR = ximisc.getDRfactors(fbase)
    else:
      DRfac = DRfacinfo[0]
      fixRR = DRfacinfo[1]


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
  dlogs = (np.log(s1d[1:]/s1d[:-1])).mean()
  if (np.fabs(mu1d[1:]-mu1d[:-1] - dmutmp) > 0.0001*dmutmp).any():
    ## not linear binning!
    ## this shouldn't be the case for a s,mu grid.
    assert 0==1
    assert np.fabs(dmutmp - dmu) < 1.0e-4

  if (np.fabs(s1d[1:]-s1d[:-1] - ds) < 0.0001*ds).all():
    logsopt = 0
    rglowedge = rg.copy()
    rglowedge = rglowedge-0.5*ds
    ## added specifically for xiopt
    rghighedge = rg.copy()
    rghighedge = rghighedge+0.5*ds
  if (np.fabs(np.log(s1d[1:])-np.log(s1d[:-1]) - dlogs) < 0.0001*dlogs).all():
    logsopt = 1
    rglowedge = rg.copy()
    rglowedge = rglowedge*np.exp(-0.5*dlogs)
    ## added specifically for xiopt
    rghighedge = rg.copy()
    rghighedge = rghighedge*np.exp(0.5*dlogs)

  assert logsopt == 0 or logsopt == 1
  if xiopt == 1:
    Vg = rghighedge**3 - rglowedge**3  ## just need something proportional to volume.
    Vzerochk = Vg.min()*0.1

  ## cut at edge of bin.
  xx = np.where(rglowedge*(1-(mug+0.5*dmu)**2)**0.5 < rperpcut)[0]
  if(rperpcut < 0.):  assert len(xx) == 0
  if len(xx) > 0:
    if periodicopt == 1:
      Vg[xx] = 0.
      Ng[xx] = 0.
    elif xiopt == 1:
      Vg[xx] = 0.
      xig[xx] = 0.

    else:
      DDg[xx] = 0.
      DRg[xx] = 0.
      RRg[xx] = 0.

  mymask = np.zeros(len(rg),dtype='int')
  mymask[xx] = 1

  ## tmp!  print a mask file.
  if(printmaskopt == 1):
    print 'yoyoyo opening masktmp.dat'
    ofpmask = open('masktmp.dat','w')
    for i in range(len(mymask)):
      ofpmask.write('%d\n' % (mymask[i]))
    ofpmask.close()

  if(printmaskopt == 2):
#### nevermind, let's print below the boundaries of each s,mu bin.
#    print 'yoyoyo opening masktmp2.dat'
#    ofpmask = open('masktmp2.dat','w')
    mumaxlist = np.zeros(len(s1d)) + 1.
    for qi in range(len(s1d)):
      rgval = s1d[qi]
      qq = np.where((rg == rgval) & (mymask == 1))[0]
      ww = np.where((rg == rgval))[0]
      assert len(ww) > 0
      if(len(qq) == 0):
        pass
#        ofpmask.write('%e %e %e\n' % (rgval,rglowedge[ww[0]],1.0))
      else:
#        ofpmask.write('%e %e %e\n' % (rgval,rglowedge[ww[0]],mug[qq].max()))
        mumaxlist[qi] = mug[qq].min()-0.5*dmu
### testing.
#      print s1d[qi], mumaxlist[qi]
      #ofpmask.write('%d\n' % (mymask[i]))
#    ofpmask.close()
    
  
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

  if printmaskopt == 2:
    ofpmask = open('masktmp2.dat','w')

  ## check mask if I'm going to cut more.
  ristart = 0
  badlist = np.zeros(nrcut,dtype='int')
  for i in range(nrcut):
    if binfile is not None:
      dfacs = rjoin[i]
    if (mymask[ristart*nmu:(ristart+dfacs)*nmu] == 1).all():
      badlist[i] = 1
    ristart = ristart + dfacs

  nrcutkeep = nrcut-len(np.where(badlist == 1)[0])

  xi = np.zeros([nell,nrcutkeep],dtype=float)
  svec = np.zeros([nell,nrcutkeep],dtype=float)

  rcen = np.zeros(nrcut,dtype=float)
  ristart = 0
  xiindx = 0
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
      if periodicopt == 1:
        if(ishort == 0):
          mymu = ximisc.downsample1d(mug[i1:i2],dfacmu)
          myVg = ximisc.downsample1dsum(Vg[i1:i2],dfacmu)
          myNg = ximisc.downsample1dsum(Ng[i1:i2],dfacmu)
        else:
          myVg = myVg + ximisc.downsample1dsum(Vg[i1:i2],dfacmu)
          myNg = myNg + ximisc.downsample1dsum(Ng[i1:i2],dfacmu)
      elif xiopt == 1:
        if(ishort == 0):
          mymu = ximisc.downsample1d(mug[i1:i2],dfacmu)
          myxig = ximisc.downsample1dvweight(xig[i1:i2],Vg[i1:i2],dfacmu)  #perform volume weighted sum. 
          myVg = ximisc.downsample1dsum(Vg[i1:i2],dfacmu)
        else:
          myxig = myxig + ximisc.downsample1dvweight(xig[i1:i2],Vg[i1:i2],dfacmu) # perform volume weighted sum, divide out at the end.
          myVg = myVg + ximisc.downsample1dsum(Vg[i1:i2],dfacmu)

      else:
        if(ishort == 0):
          mymu = ximisc.downsample1d(mug[i1:i2],dfacmu)
          mydd = ximisc.downsample1dsum(DDg[i1:i2],dfacmu)
          mydr = ximisc.downsample1dsum(DRg[i1:i2],dfacmu)
          myrr = ximisc.downsample1dsum(RRg[i1:i2],dfacmu)
        else:
          mydd = mydd + ximisc.downsample1dsum(DDg[i1:i2],dfacmu)
          mydr = mydr + ximisc.downsample1dsum(DRg[i1:i2],dfacmu)
          myrr = myrr + ximisc.downsample1dsum(RRg[i1:i2],dfacmu)

    if periodicopt == 1:
      yy = np.where(myVg < Vzerochk)[0]
      xitmp = myNg/myVg-1.
      xitmp[yy] = 0.

    elif xiopt == 1:
      yy = np.where(myVg < Vzerochk)[0]
      xitmp = myxig/myVg
      xitmp[yy] = 0.
      ## for now, nothing to 0 out (?)

    else:  ##DR stuff
      zz = np.where(myrr < smallRRcut)[0]
      if(len(zz) > 0):
        print 'using smallRRcut!  here are details',i,rcen[i],smallRRcut
        myrr = myrr.mean()

      yy = np.where(myrr < 0.01)[0]
      xitmp = (mydd-mydr*DRfac)/myrr/DRfac**2/fixRR**2+1.
      xitmp[yy] = 0.  ## fix 0'd out regions

      if(len(yy) > 0 and 0==1):
        print 'ack why zerod out regions?'
        print len(yy)
        print i, dfacs, dfacmu
        print myrr
        print mymu

    for ell in np.arange(0,nell*2,2):
      if(badlist[i] == 0):
        xi[ell/2,xiindx] = (xitmp*ximisc.legendre(ell,mymu)).sum()*dmudown*(2.*ell+1.)
        svec[ell/2,xiindx] = rcen[i]

    if printmaskopt == 2:
      ofpmask.write('%e %e %e %e %e ' % (rcen[i], s1d[ristart]*np.exp(-0.5*dlogs), s1d[riend]*np.exp(0.5*dlogs), mumaxlist[ristart],mumaxlist[riend]))
      for ell in np.arange(0,nell*2,2):
        if badlist[i] == 0:
          ofpmask.write('%e ' % (xi[ell/2,xiindx]))
        else:
          ofpmask.write('%e ' % (0.0))

      ofpmask.write('\n')

    ristart = ristart + dfacs

    if(badlist[i] == 0):
      xiindx += 1

  return xiell(xiellfname=fDR,sxilist=[svec, xi],icovfname=icovfname)

## stealing code from publicly released code
## http://www.sdss3.org/science/boss_publications.php
## also on howdiedoo ~/dr9reidetalpublic

def invcovcombine(xia,xib,outfilebase,xi0istart=None,xi0iend=None,xi2istart=None,xi2iend=None,unbiasicovfac=1.,unbiascovfac=1.):
  """
  take two correlation functions and combine them with inverse cov weighting,
  using only a subset of bins if desired.
  set istart,iend = -1 if you want to exclude that multipole
  first bin is indexed by 0
  to include all bins, use: combine.py 0 len(xia.xi0)-1 0 len(xia.xi2)-1
  or leave default None (will set to full length).
  cov and xi will be written using outfilebase, then read back in to xiell to return an object.
  unbias factors not written to the files, just passed along to the xiell.xiell call.
  """

  print 'this function was written on july 29 and kind of tested against.'
  print 'see /home/howdiedoo/boss/cosmoxi2dfitalphadr11/data/ladodr11v1/combinestests'
  print 'you should still use it with caution.'

  if (len(xia.xi0) != len(xib.xi0)):
    print 's0 binnings dont agree, returning none!'
    return None
  if (np.fabs(xia.s0 - xib.s0) > 2.0e-6).any():
    print 's0 binnings dont agree, returning none!'
    return None
  if (xi2istart is not None) and xi2istart >= 0: ## that means you want to use xi2
    if (np.fabs(xia.s2 - xib.s2) > 2.0e-6).any():
      print 's2 binnings dont agree, returning none!'
      return None

  ## everything is coded currently so that xi0 and xi2 binning is the same, so test for this.
  if (xi2istart is not None) and xi2istart >= 0: ## that means you want to use xi2
    ellmax = 2
    if len(xia.s0) != len(xia.s2):
      print 's0 and s2 are different, returning none!'
      return None
    if (np.fabs(xia.s0 - xia.s2) > 2.0e-6).all():
      print 's0 and s2 are different, returning none!'
      return None
  else:
    ellmax = 0    

  imin = 0
  imax = len(xia.xi0)-1
  nbinstot = imax - imin + 1
  assert nbinstot == len(xia.xi0)
  assert nbinstot == len(xia.xi2)

  if(xi0iend == -1):
    assert xi0istart == -1
    nxi0 = 0
  else:
    nxi0 = xi0iend-xi0istart+1
  if((xi2iend is None) or (xi2iend == -1)):
    assert xi2istart == -1
    nxi2 = 0
  else:
    nxi2 = xi2iend-xi2istart+1

  ncuttot = nxi0 + nxi2

  mynell = 1
  if nxi2 > 0:
    mynell = 2

  if(xi0istart == -1):
    xi0iend = -1
  if(xi2istart == -1):
    xi2iend = -5*nbinstot

  if xi0istart == -1 or xi0istart is None:
    print 'xi2 only not currently supported, need to rewrite xiell.xiell() for that!'
    return None

  ncov1d = len(xia.cov[:,0])
  if ncov1d != nbinstot*2 and nxi2 > 0:
    print 'cov make no sense, missing xi2 in cov?',ncov1d,nbinstot,nbinstot*2
    return None
  if ncov1d != nbinstot*2 and ncov1d != nbinstot:
    print 'cov make no sense, bins not aligned?',ncov1d,nbinstot,nbinstot*2
    return None

  ## do several cases. 
  ## case 1: cov is xi0 only

  ## this part copied from combine.py
  
  scut = np.zeros(ncuttot)
  xiacut = np.zeros(ncuttot)
  xibcut = np.zeros(ncuttot)
  covacut = np.zeros([ncuttot,ncuttot])
  covbcut = np.zeros([ncuttot,ncuttot])
  
  ## picks out the relevant rows/columns from xi and cov
  i1=0
  for j1 in range(ncov1d):
    if(not (((j1 >= xi0istart) and (j1 <= xi0iend)) or ((j1-nbinstot >= xi2istart) and (j1-nbinstot <= xi2iend)))):
      continue
    if(j1 < nbinstot):
      scut[i1] = xia.s0[j1]
      xiacut[i1] = xia.xi0[j1]
      xibcut[i1] = xib.xi0[j1]
    else:
      scut[i1] = xia.s2[j1-nbinstot]
      xiacut[i1] = xia.xi2[j1-nbinstot]
      assert xia.xilong[j1] == xia.xi2[j1-nbinstot]
      xibcut[i1] = xib.xi2[j1-nbinstot]
      assert xib.xilong[j1] == xib.xi2[j1-nbinstot]
    i2=0
    for j2 in range(ncov1d):
      if(not (((j2 >= xi0istart) and (j2 <= xi0iend)) or ((j2-nbinstot >= xi2istart) and (j2-nbinstot <= xi2iend)))):
        continue
      covacut[i1,i2] = xia.cov[j1,j2]
      covbcut[i1,i2] = xib.cov[j1,j2]
      i2 += 1
    i1 += 1

  icovacut = (np.matrix(covacut)).I
  icovbcut = (np.matrix(covbcut)).I
  icovtot = icovacut + icovbcut
  covtot = (icovtot).I
  xitotcut = covtot*(icovacut*np.transpose(np.matrix(xiacut)) + icovbcut*np.transpose(np.matrix(xibcut)))

  xitotcut = np.array(xitotcut[:,0])

  cov.printcov(np.array(covtot),outfilebase+'.cov')
  cov.printcov(np.array(icovtot),outfilebase+'.icov')
  icovfname = outfilebase+'.icov'
  ## write the cov to a file in the usual format, and read it back in so that everything is exactly reproducible.
  xinew = xiell(sxilist=[scut.reshape(mynell,len(scut)/mynell),xitotcut.reshape(mynell,len(xitotcut)/mynell)],icovfname=icovfname,unbiasicovfac=unbiasicovfac, unbiascovfac=unbiascovfac)

  xinew.printxielllong(outfilebase+'.xi')

  return xinew 
   


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



