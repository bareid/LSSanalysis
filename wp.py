import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

class wp:
  def __init__(self,wpfname,icovfname=None):
    try:
      rsig, wp = np.loadtxt(wpfname,unpack=True,usecols=[0,1])
      self.fname = [wpfname]
      self.nrsig = len(rsig)
      self.rsig = rsig
      self.wp = wp
      ## comput binning.
      self.logopt = -1 ## non-simple binning.
      drsig = (self.rsig[1:]-self.rsig[:-1]).mean()
      if (np.fabs(self.rsig[1:]-self.rsig[:-1] - drsig) < 0.0001*drsig).all():
        self.logopt = 0
        self.drsig = drsig
      dlogrsig = (np.log(self.rsig[1:]/self.rsig[:-1])).mean()
      if (np.fabs(np.log(self.rsig[1:])-np.log(self.rsig[:-1]) - dlogrsig) < 0.0001*dlogrsig).all():
        self.logopt = 1
        self.dlogrsig = dlogrsig
      
    except:
      print 'bad wp file.'
      self = None

    if(self is not None):
      ## let's finde out if this is a data of theory curve.
      fbase = wpfname.split('.'+wpfname.split('.')[-1])[0]
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
        assert len(icov[:,0]) == self.nrsig
## compute diag errors.
        cov = icov.I
        diagerr = np.array([np.sqrt(cov[i,i]) for i in range(self.nrsig)])
        self.diagerr = diagerr 

  def __str__(self):
    mystr = "%d sized wp from [%s]. Range = [%f %f], logopt = %d\n" % (self.nrsig,self.fname,self.rsig.min(),self.rsig.max(),self.logopt)
    return mystr

  def checkbinning(self,other):
    """
    Returns 0 if binning is the same, 1 otherwise
    """
    if(self.nrsig != other.nrsig):
      return 1
    if((np.fabs(self.rsig - other.rsig) > 2.0e-4).any()):
      return 1
    return 0


  def __add__(self,other):
    if(self.checkbinning(other) == 0):
      wpnew = copy.deepcopy(self)
      wpnew.wp = (self.wp*self.weight + other.wp*other.weight)/(self.weight+other.weight)
      wpnew.weight = (self.weight+other.weight)
      ## this needs to be generalized if you want to add ones that have already been added.
      wpnew.fname = self.fname + other.fname
    else:
      print 'add fail!  binning misaligned'
      wpnew = None
    return wpnew

  def wpinterp(self,rsigarr):
    ## do interpolation in log space.
    if ((type(rsigarr) is type(0.0)) | (type(rsigarr) is type(0))):
      rsigarr = np.array([rsigarr])

    ## make sure it's type numpy
    rsigarr = np.array(rsigarr)

    if(self.logopt == -1):
## need to write this piece later.
      return 1.

    if(self.logopt == 0):
      ix = (rsigarr-self.rsig.min())/self.drsig 
    if(self.logopt == 1):
      ix = (np.log(rsigarr/self.rsig.min()))/self.dlogrsig 
    tmp = np.where(ix < 0)[0]
    ix[tmp] = 0
    tmp = np.where(ix > (self.nrsig)-1)[0]
    ix[tmp] = self.nrsig-1.001
    iix = np.array([int(val) for val in ix])
    assert (iix >= 0).all()
    assert (iix <= self.nrsig-2).all()
    fx = ix -iix
    assert (fx >= 0.).all()
    assert (fx <= 1.).all()


    if(self.logopt == 0):
      return self.wp[iix]*(1.-fx) + self.wp[iix+1]*fx
    if(self.logopt == 1):
      return np.exp(np.log(self.wp[iix])*(1.-fx) + np.log(self.wp[iix+1])*fx)

  def addcurve(self,ax,color='k',rppow=0,fmt=None,lbl=None):
    """
    adds a curve on an rp-wp plot.
    """
    if(self.DorT == 0): ## data, so there are errors.
      ax.errorbar(self.rsig,self.wp*self.rsig**rppow,yerr=self.diagerr*self.rsig**rppow,color=color,ecolor=color,fmt=fmt,label=lbl)
    else:
      ax.plot(self.rsig,self.wp*self.rsig**rppow,label=lbl,color=color)

  def makeplot(self,ax=None,color='k',span=None,logxopt=1,logyopt=1,rppow=0,fmt=None,lbl=None,customxticks=[0.5,1.0,5.0,10.,20.,30.]):
    """
    Makes a standard wp one-d plot.
    log[x-y]opt specify if you want log or linear binning on the axis.
    rppow means multiply y coordinate by rsig**rppow.
    """

    if ax is None:
      ff = plt.figure(figsize=[6,6])
      ax=ff.add_subplot(1,1,1)
    else:
      ff = None

    self.addcurve(ax,color=color,rppow=rppow,lbl=lbl,fmt=fmt)

    if span is None:
      span = [self.rsig.min()*0.9, self.rsig.max()*1.1, (self.wp*self.rsig**rppow).min()*0.9, (self.wp*self.rsig**rppow).max()*1.1]
    if(logxopt == 1):
      ax.set_xscale('log')
    if(logyopt == 1):
      ax.set_yscale('log')


    ax.axis(span)
    ax.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=20)
    if(np.fabs(rppow) > 0.01):
      ax.set_ylabel(r'$r_{\sigma}^{%.1f} w_p(r_{\sigma})$' % (rppow),fontsize=16)
    else:
      ax.set_ylabel(r'$w_p(r_{\sigma}) \, [h^{-1} {\rm Mpc}]$',fontsize=20)

    ## set up some custom x axis labels if desired.
    if(customxticks is not None):
      ax.xaxis.set_major_locator(plt.FixedLocator(customxticks))
      ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
      ax.xaxis.set_ticks_position('bottom')
      ax.tick_params(axis='x',reset=False,which='both',length=8,width=2)

    return ff, ax

  def chi2(self,other):
    """
    compute chi2 of self with other.  self is the one with the inverse cov.
    """
    diff = np.matrix(self.wp-other.wp)
    if self.DorT == 0 and other.DorT == 1:
      return ((diff*self.icov) * (diff.T))[0,0]
    elif self.DorT == 1 and other.DorT == 0:
      return ((diff*other.icov) * (diff.T))[0,0]
    else:
      return 0. ## can't compute chi2 without an inverse covariance matrix.



if __name__ == "__main__":

  wp = wp(wpfname="../boss/bethalexie/data/wpNSmean-dr10_fbworkCHK-nofkp-nocp-smallbins-diagerrs.dat.corr",icovfname="../boss/bethalexie/data/icov-dr10_fbworkCHK-nofkp-nocp-smallbins.dat")
  print wp
  print wp.rsig
  print wp.wp
  
  print wp.nrsig
  print wp.logopt
  if(wp.logopt == 1): print wp.dlogrsig 
  if(wp.logopt == 0): print wp.drsig 
  if(wp.logopt == -1): print 'weird binning'


