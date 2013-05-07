import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

class xiell:
  def __init__(self,xiellfname,icovfname=None,ellmax=-1):
    """
    If you only want to read in some of the ell values, use ellmax.
    Otherwise reads in all the multipoles in the file.
    """
#    try:
    if(0==0):
      self.fname = [xiellfname]
      ## accept 2 formats.
## later extend format to include xi0,2,4?
      aa = np.loadtxt(xiellfname)
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
      self.xilong = aa[:,1]
      self.ndata = len(aa[:,1])


      if(nxi > 1):
        nell = nxi
        if(ellmax != -1 and (ellmax == 0 or ellmax == 2 or ellmax == 4)):
          nell = ellmax/2+1
        
        xi = np.zeros([nell,ns])
        svec = np.zeros([nell,ns])
        xilong = []
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

#    except:
    else:
      print 'bad xiell file',xiellfname
      self = None

    if(self is not None):
      ## let's finde out if this is a data of theory curve.
      fbase = xiellfname.split('.'+xiellfname.split('.')[-1])[0]
      try:
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
      self.fig = plt.figure(figsize=[6,6])
      ax=self.fig.add_subplot(1,1,1)

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

    return ax

  def makediffplot(self,other,ell=0,ax=None,color='k',span=None,logxopt=1,logyopt=0,spow=1,fmt=None,lbl=None,ecolor='k'):
    """
    Makes a standard xiell 1d plot, but plot difference between data (self) and theory (other).
    log[x-y]opt specify if you want log or linear binning on the axis.
    spow means multiply y coordinate by s**spow.
    """

    if ax is None:
      self.fig = plt.figure(figsize=[6,6])
      ax=self.fig.add_subplot(1,1,1)

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

    return ax




if __name__ == '__main__':
  print 'hi'






