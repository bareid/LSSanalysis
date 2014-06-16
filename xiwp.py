import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ximisc
import cov
import xiell
import wp

class xiwp:
  def __init__(self,xiellin=None,wpin=None,xiellwpfname=None,icovfname=None):
    if xiellin is not None:
      self.xiell = xiellin
    if wpin is not None:
      self.wp = wpin

    if xiellwpfname is not None:
      self.xiellwpfname = xiellwpfname

    if xiellin is not None and wpin is not None:
      self.ntot = self.wp.nrsig + self.xiell.ndata
      self.xiwp = np.concatenate((self.xiell.xilong, self.wp.wp))

      if len(self.xiwp) != self.ntot:
        print 'vector length mismatch!'
        self = None

    else:
      if xiellwpfname is None:
        self = None
        return

      try:
        rx, self.xiwp = np.loadtxt(xiellwpfname,usecols=[0,1],unpack=True)
        self.ntot = len(self.xiwp)
      except:
        self = None

      if self is not None:
        try:
          ## find xiell/wp split, create 
          xx = np.where(rx[:-1] > rx[1:])[0]
          ## i *think* this should work, double check if you choose some wacko binning!
          xwp = xx[-1]+1
          nell = len(xx)
          svec = rx[:xwp].reshape(nell,len(rx[:xwp])/nell)
          xivec = self.xiwp[:xwp].reshape(nell,len(rx[:xwp])/nell)
          self.xiell = xiell.xiell(sxilist=[svec,xivec])
          self.wp = wp.wp(rpwplist=[rx[xwp:], self.xiwp[xwp:]])

        except: # oh well.
          pass

    ## tmp.
#    self.xiell = xiell.xiell(sxilist=[rx[:xwp],self.xiwp[:xwp]])
#    self.wp = wp.wp(rpwplist=[rx[xwp:], self.xiwp[xwp:]])

  
    if self is not None:
      if icovfname is None:
        self.DorT = 1
      else:
        self.DorT = 0
        icov = np.matrix(np.loadtxt(icovfname))
        self.icov = icov
        try:
          cov = icov.I
          diagerr = np.array([np.sqrt(cov[i,i]) for i in range(self.ntot)])
          self.diagerr = diagerr
        except:
          diagerr = np.zeros(self.ntot) + 1.
          self.diagerr = diagerr
  
        if not len(icov[:,0]) == self.ntot:
          self = None

  def __str__(self):
    mystr = "%d sized xiwp." % self.ntot

  def printxiwp(self,outfname):
    ofp = open(outfname,'w')
    if self.DorT == 1:
      for i in range(self.ntot):
        ofp.write('%e\n' % (self.xiwp[i]))
    else:
      for i in range(self.ntot):
        ofp.write('%e %e\n' % (self.xiwp[i], self.diagerr[i]))
    ofp.close()

  def chi2(self,other):
    """
    compute chi2 of self with other.
    """
    diff = np.matrix(self.xiwp-other.xiwp)
    if self.DorT == 0 and other.DorT == 1:
      return ((diff*self.icov) * (diff.T))[0,0]
    elif self.DorT == 1 and other.DorT == 0:
      return ((diff*other.icov) * (diff.T))[0,0]
    elif self.DorT == 0 and other.DorT == 0:
      return ((diff*self.icov) * (diff.T))[0,0]
    else:
      return 0. ## can't compute chi2 without an inverse covariance matrix.

  def chi2info(self,other):
    """
    For non-diagonal cov, this helps show you which differences between model and theory contribute the most
    to chi2.
    """
    if isinstance(other,xiwp):
      diff = np.matrix(self.xiwp-other.xiwp)
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
        diff = np.matrix(self.xiwp - other)
        myicov = self.icov
      except:
        return 0.
    mysum = 0.
    for i in range(len(self.xiwp)):
      print i,diff[0,i]*((diff*myicov)[0,i])
      mysum += diff[0,i]*((diff*myicov)[0,i])
    assert np.fabs(mysum - self.chi2(other)) < 0.001
    return mysum



  def makeplot(self):
    """
    Coming soon!  Actually you should use makexiwpplot which inputs wp and xiell (rather than a xiMwp object).
    """


  def sanitycheckicov(xiellwcov,wpwcov):
    """
    Input wp and xiell with separate icovs derived in the same manner as self.
    we'll make sure everything is in correct order by checking agreement for 
    diagonal elts of covs.
    """
    for i in range(self.xiell.ndata):
      assert np.fabs(self.diagerr[i] - self.xiell.diagerr[i])/self.diagerr[i] < 1.0e-4
    for i in range(self.wp.nrsig):
      assert np.fabs(self.diagerr[self.xiell.ndata+i] - self.wp.diagerr[i])/self.diagerr[self.xiell.ndata+i] < 1.0e-4
    print 'passed sanity check on diagerrs for icov!!'
    return 0

def makeplotxiwpplot(xwp,xxiell,xxiell2=None,axlist=None,colorlist=None,colorlist2=None,rppow=0,wpfmt=None,wplbl=None,wpcustomxticks=[0.5,1.0,5.0,10.,20.,30.],spow=1,xifmt=None,xilbl=None,xilbl2=None,diffopt=0,xifmt2=None):

  if axlist is None:
    if xxiell2 is not None:
      f = plt.figure(figsize=[17,5.5])
      axlist = []
      axlist.append(f.add_subplot(131))
      axlist.append(f.add_subplot(132))
      axlist.append(f.add_subplot(133))
    else:
      f = plt.figure(figsize=[12,5.5])
      axlist = []
      axlist.append(f.add_subplot(121))
      axlist.append(f.add_subplot(122))

  else:
    f = None

  if diffopt == 1 and xxiell2 is None:
    diffopt = 0

  if colorlist is None:
    colorlist = ['k','k','k']
  if colorlist2 is None:
    colorlist2 = [colorlist[1], colorlist[2]]
  if xifmt2 is None:
    xifmt2 = xifmt
  if xilbl2 is None:
    xilbl2 = xilbl

  xwp.makeplot(ax=axlist[0],color=colorlist[0],rppow=rppow,fmt=wpfmt,lbl=wplbl,customxticks=wpcustomxticks)
  if diffopt == 1:
    xxiell.makediffplot(xxiell2,ax=axlist[1],color=colorlist[1],ell=0,spow=spow,fmt=xifmt,lbl=xilbl)
    xxiell.makediffplot(xxiell2,ax=axlist[2],color=colorlist2[0],ell=2,spow=spow,fmt=xifmt2,lbl=xilbl2)
  
  else:
    print 'acko',xifmt2
    xxiell.makeplot(ax=axlist[1],color=colorlist[1],ell=0,spow=spow,fmt=xifmt,lbl=xilbl)
    xxiell.makeplot(ax=axlist[1],color=colorlist2[0],ell=2,spow=spow,fmt=xifmt2,lbl=xilbl2)
    if xxiell2 is not None and len(axlist) > 2 and len(colorlist) > 2 and len(colorlist2) > 1:
      xxiell2.makeplot(ax=axlist[2],color=colorlist[2],ell=0,spow=spow,fmt=xifmt,lbl=xilbl)
      xxiell2.makeplot(ax=axlist[2],color=colorlist2[1],ell=2,spow=spow,fmt=xifmt2,lbl=xilbl2)

  return f, axlist 

if __name__ == "__main__":
  pass


