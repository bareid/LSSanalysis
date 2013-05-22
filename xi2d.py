#!/usr/bin/python
# Filename: xi2d.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

def getflip(xi,n1d):
  xiflip = xi.copy()
  for i in range(n1d):
    for j in range(n1d):
      xiflip[i,j] = xi[j,i]
  return xiflip

class xi2d:
  def __init__(self,xi2dfname):
    try:
      rsig, rpi, xi = np.loadtxt(xi2dfname,unpack=True)
      n1d = int(np.sqrt(len(rsig)))
      if(np.fabs(n1d**2-len(rsig)) > 0):
        print 'Have not implemented non-square xi2ds. Try back later!'
        sys.exit(1)

      ## make it a list, could have come from more than one place.
      self.fname = [xi2dfname]
      self.n1d = n1d
      self.nrsig = n1d
      self.nrpi = n1d
      self.ntot = len(rsig)
      self.rsig = rsig
      self.rpi = rpi
      self.xi = xi
      self.rsig1d = rsig.reshape(n1d,n1d)[:,0]
      self.rpi1d = rpi.reshape(n1d,n1d)[0,:]

      ## compute linear spacings.
      drsig = (self.rsig1d[1:]-self.rsig1d[:-1]).mean()
      drpi = (self.rpi1d[1:]-self.rpi1d[:-1]).mean()
      if (np.fabs(self.rsig1d[1:]-self.rsig1d[:-1] - drsig) > 0.0001*drsig).any():
        ## not linear binning!
        self.drsig = 0.
      else:
        self.drsig = drsig 
      if (np.fabs(self.rpi1d[1:]-self.rpi1d[:-1] - drpi) > 0.0001*drpi).any():
        ## not linear binning!
        self.drpi = 0.
      else:
        self.drpi = drpi

    except:
#    if(0==1):
      print 'bad xi2d file.'
      self = None

    if(self):
      ## let's finde out if this is a data of theory curve.
      fbase = xi2dfname.split('.'+xi2dfname.split('.')[-1])[0]
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
#        print 'assuming this is a theory xi2d rather than data'
        self.weight = 1.0



  def __str__(self):
    mystr = "(%d, %d) sized xi2d from [%s].  Range = [[%f, %f], [%f, %f]]\n" % (self.nrsig,self.nrpi,", ".join(self.fname),self.rsig.min(), self.rsig.max(),self.rpi.min(),self.rpi.max())
    return mystr

  def checkbinning(self,other):
    """
    Returns 0 if binning is the same, 1 otherwise
    """
    if(self.n1d != other.n1d):
      return 1
    if(self.ntot != other.ntot):
      return 1
    if(self.nrsig != other.nrsig):
      return 1
    if(self.nrpi != other.nrpi):
      return 1
    if((np.fabs(self.rsig - other.rsig) > 2.0e-4).any()):
      return 1
    if((np.fabs(self.rpi - other.rpi) > 2.0e-4).any()):
      return 1
    return 0

  def __add__(self,other):
    if(self.checkbinning(other) == 0):
      xinew = copy.deepcopy(self)
      xinew.xi = (self.xi*self.weight + other.xi*other.weight)/(self.weight+other.weight)
      xinew.weight = (self.weight+other.weight)
      ## this needs to be generalized if you want to add ones that have already been added.
      xinew.fname = self.fname + other.fname
    else:
      print 'add fail!  binning misaligned'
      xinew = None
    return xinew

  def xiinterp(self, rsigarr,rpiarr):
    """ Quick n dirty bilinear interpolation in xi
        Make this operate on input array of values.
        should work for a single float value as well.
    """

    if ((type(rsigarr) is type(0.0)) | (type(rsigarr) is type(0))):
      rsigarr = np.array([rsigarr])
      rpiarr = np.array([rpiarr])

    ## make sure it's type numpy
    rsigarr = np.array(rsigarr)
    rpiarr = np.array(rpiarr)

    ix = (rsigarr-self.rsig1d.min())/self.drsig
    tmp = np.where(ix < 0)[0]
    ix[tmp] = 0
    tmp = np.where(ix > (self.nrsig)-1)[0]
    ix[tmp] = self.nrsig-1.001

    iy = (rpiarr-self.rpi1d.min())/self.drpi
    tmp = np.where(iy < 0)[0]
    iy[tmp] = 0
    tmp = np.where(iy > (self.nrpi)-1)[0]
    iy[tmp] = self.nrpi-1.001

    iix = np.array([int(val) for val in ix])
    iiy = np.array([int(val) for val in iy])
    assert (iix >= 0).all()
    assert (iix <= self.nrsig-2).all()
    assert (iiy >= 0).all()
    assert (iiy <= self.nrpi-2).all()
    xi11 = self.xi[iix*(self.nrpi) + iiy]
    xi21 = self.xi[(iix+1)*(self.nrpi) + iiy]
    xi12 = self.xi[iix*(self.nrpi) + iiy+1]
    xi22 = self.xi[(iix+1)*(self.nrpi) + iiy+1]

    ## fraction of a bin left over.
    fx = ix -iix
    fy = iy -iiy
    assert (fx >= 0.).all()
    assert (fx <= 1.).all()
    assert (fy >= 0.).all()
    assert (fy <= 1.).all()

    xival = xi11*(1.-fx)*(1.-fy) + \
            xi21*fx*(1.-fy) + \
            xi12*(1-fx)*(fy) + \
            xi22*fx*(fy)

    if(len(xival) == 1):
      xival = xival[0]
    return xival

  def getclevels(self,ncontours=None,cratio=0.7,dc=0.02,ximin=1.0e-3,ximax=8.0,logopt=1,useself=False):
    """
    Function to intelligently guess contour level values.
    If useself = False, uses default or input values.
    If True, uses min/max of self.xi
    """

    if(useself == True):
      if(logopt == 1):
        ximax = self.xi.max()
        xx = np.where(xi > 0.)
        ximin = (self.xi[xx]).min()
        if(ximax < 0.): logopt = 0
      if(logopt != 1):
        ximax = self.xi.max()
        ximin = self.xi.min()
        

    if ncontours is not None and ncontours > 2:
      if(logopt == 1):
        cratio = np.exp(np.log(ximin/ximax)/float(ncontours-1))
        dc = 0
      else:
        cratio = 1.0
        dc = (ximax-ximin)/float(ncontours-1)

    cval = ximax
    tclevlist = [cval]
    if(logopt == 1):
      ximin = np.fabs(ximin)
      cnt = 0
      while(cval > ximin*0.999 and cnt < 1000):
        cval = cval*cratio
        tclevlist.append(cval)
        cnt += 1
    else:
      cnt = 0
      while(cval > ximin-N.fabs(0.001*ximin) and cnt < 1000):
        cval = cval-dc
        tclevlist.append(cval)
        cnt += 1
    if ncontours is not None:
      if len(tclevlist) > ncontours:
        tclevlist = tclevlist[:ncontours]
    return tclevlist

  def symmetrize(self,symmetrizeopt=1):
    n1d = self.n1d
    if(symmetrizeopt == 1):
      n1dS = 2*self.n1d
      off = self.n1d
      rsigS = np.zeros([n1dS,n1dS])
      rpiS = np.zeros([n1dS,n1dS])
      xiS = np.zeros([n1dS,n1dS])
      for i in range(0,self.n1d):
        for j in range(0,self.n1d):
          rsigS[off+i, off+j] = self.rsig.reshape(n1d,n1d)[i,j]
          rsigS[off-i-1, off+j] = -self.rsig.reshape(n1d,n1d)[i,j]
          rsigS[off-i-1, off-j-1] = -self.rsig.reshape(n1d,n1d)[i,j]
          rsigS[off+i, off-j-1] = self.rsig.reshape(n1d,n1d)[i,j]

          rpiS[off+i, off+j] = self.rpi.reshape(n1d,n1d)[i,j]
          rpiS[off-i-1, off+j] = self.rpi.reshape(n1d,n1d)[i,j]
          rpiS[off-i-1, off-j-1] = -self.rpi.reshape(n1d,n1d)[i,j]
          rpiS[off+i, off-j-1] = -self.rpi.reshape(n1d,n1d)[i,j]

          xiS[off+i, off+j] = self.xi.reshape(n1d,n1d)[i,j]
          xiS[off-i-1, off+j] = self.xi.reshape(n1d,n1d)[i,j]
          xiS[off-i-1, off-j-1] = self.xi.reshape(n1d,n1d)[i,j]
          xiS[off+i, off-j-1] = self.xi.reshape(n1d,n1d)[i,j]
      rsigS1d = rsigS[:,0]
      rpiS1d = rpiS[0,:]

    else:
      rsigS = self.rsig
      rpiS = self.rpi
      xiS = self.xi
      n1dS = self.n1d
      rsigS1d = self.rsig1d
      rpiS1d = self.rpi1d
  

    return rsigS1d, rpiS1d, xiS, n1dS
  
  def addcontour(self,ax,symmetrizeopt=1,clevlist=[],color='k'):
    """
    adds a contour plot of the current object to axis ax.
    """
    if(clevlist == []):
      clevlist=self.getclevels()

    rsigS1d, rpiS1d, xiS, n1dS = self.symmetrize(symmetrizeopt)
    xiflip = getflip(xiS.reshape(n1dS,n1dS),n1dS)
    cc = ax.contour(rsigS1d,rpiS1d,xiflip,clevlist,colors=color,linewidths=2)
    return cc 

  ## things look screwy if you don't symmetrize, so we don't allow unsymmetrized for now.
  ## remap options do some sort of transformation on xi to make the contours nicer.
  ## so far 
  def adddensity(self,ax,remapopt=1):
    """
    Make a density plot.
    remapopt = 0 does nothing to xi
    remapopt = 1 makes density plot of log(xi) [this will die if xi < 0!]
    add more options as needed.
    """
    rsigS1d, rpiS1d, xiS, n1dS = self.symmetrize(1)
    if(remapopt == 1):
      xiremap = (np.log10(xiS)/np.log10(self.xi.max()/self.xi.min()))
    else:
      xiremap = xiS
    xiflip = getflip(xiremap.reshape(n1dS,n1dS),n1dS)
    xl=self.rsig.max()
    yl=self.rpi.max()
    im = ax.imshow(xiflip,extent=[-xl,xl,-yl,yl])
    return im 


  def makecontourplot(self,ax=None,symmetrizeopt=1,clevlist=[],color='k',span=None):
    """
    Make a contour plot owned by this object if ax=None, returns the axis.
    Otherwise, just add the contour plot to input axis ax and beautify.
    Optional span = [xmin, xmax, ymin, ymax]
    """
    if ax is None:
## we want aspect ratio to be 1.!!
      ff = plt.figure(figsize=[6,6])
      ax=ff.add_subplot(1,1,1)
    else:
      ff = None

    self.addcontour(ax,symmetrizeopt,clevlist,color)
    if span is None:
      if(symmetrizeopt == 1):
        span = [-self.rsig.max(), self.rsig.max(), -self.rpi.max(), self.rpi.max()]
      else:
        span = [self.rsig.min(), self.rsig.max(), self.rpi.min(), self.rpi.max()]

    ax.axis(span)
    if ff is not None:
      ax.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
      ax.set_ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
    return ff, ax
    
  def makedensityplot(self,ax=None,span=None,remapopt=1):
    """
    Make a density plot owned by this object if ax=None, returns the axis.
    Otherwise, just add the density plot to input axis ax and beautify.
    Optional span = [xmin, xmax, ymin, ymax]
    remapopt = 0 does nothing, remapopt = 1 takes a log of xi.
    """
    if ax is None:
## we want aspect ratio to be 1.!!
      ff = plt.figure(figsize=[6,6])
      ax=ff.add_subplot(1,1,1)

    self.adddensity(ax,remapopt)
    if span is None:
      span = [-self.rsig.max(), self.rsig.max(), -self.rpi.max(), self.rpi.max()]

    ax.axis(span)
    ax.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
    ax.set_ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
    return ff, ax

if __name__ == "__main__":

  if(0==1):
    myxireal = xi2d(xi2dfname="makeHJcatalogv2zspace2.cat.zspace-1.xiell.nbins900.butterfly")
    myxi = xi2d(xi2dfname="makeHJcatalogv2zspace2.cat.zspace2.xiell.nbins900.butterfly")
    myxifine = xi2d(xi2dfname="makeHJcatalogv2zspace2.cat.zspace2.xiell.nbins3600.butterflyfine")
    myxicoarse = xi2d(xi2dfname="makeHJcatalogv2zspace2.cat.zspace2.xiell.nbins225.butterflycoarse")
    plt.figure()
  #  myxireal.plotcontour(color='k')
    myxi.plotcontour(color='b')
    myxifine.plotcontour(color='g')
    myxicoarse.plotcontour(color='r')
    plt.savefig("contoursv0.png")

  ## try reading in data.
  xiN = xi2d(xi2dfname='../boss/zdistvXlogbinsompcleverLSbutterfly/outputmksamplelatestdr10v7/collidedBR-collate-cmass-dr10v7-N-FBBRang_xigrid.butterfly')
  xiS = xi2d(xi2dfname='../boss/zdistvXlogbinsompcleverLSbutterfly/outputmksamplelatestdr10v7/collidedBR-collate-cmass-dr10v7-S-FBBRang_xigrid.butterfly')
  if(xiN is not None and xiS is not None):
    print xiN.weight, xiS.weight, xiN.weight/(xiN.weight+xiS.weight)
  else:
    print 'blerg errors'
    sys.exit(1)

  xidata = xiN + xiS
  print xiN
  print xiS
  print xidata
### xi data makes sense now!
#  for i in range(xidata.ntot):
#    print xidata.rsig[i], xiN.rsig[i], xidata.rpi[i], xiN.rpi[i], xidata.xi[i], xiN.xi[i], xiS.xi[i], (xiN.xi[i]*xiN.weight+xiS.xi[i]*xiS.weight)/xidata.weight

