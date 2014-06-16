import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ximisc

class wp:
  def __init__(self,wpfname=None,icovfname=None,rpwplist=[],wpstart=-1,wpend=-1):
    """
    use wpstart,end to restrict data points to wpstart:wpend+1
    """

    if rpwplist != []:
      try:
        rsig = rpwplist[0]
        wp = rpwplist[1]
        self.fname = wpfname
      except:
        print 'bad rpwplist'
        self = None
    else: #read from a file.
      try:
        rsig, wp = np.loadtxt(wpfname,unpack=True,usecols=[0,1])
        self.fname = [wpfname]
      except:
        print 'bad wp file'
        self = None

    ## dec 30 adding wpstart, wpend
    mywpstart = wpstart
    mywpend = wpend
    if wpstart == -1:
      mywpstart = 0
    if wpend == -1:
      mywpend = len(wp)-1
    ## restrict to wpstart:wpend+1
    rsig = rsig[mywpstart:mywpend+1]
    wp = wp[mywpstart:mywpend+1]
    ## end wpstart/end

    if self is not None:
      try:
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
        print 'bad rp/wp binning or something.'
        self = None
    if(self is not None):
      ## let's finde out if this is a data of theory curve.
      try:
        fbase = wpfname.split('.'+wpfname.split('.')[-1])[0]
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
        try:
          cov = icov.I
          diagerr = np.array([np.sqrt(cov[i,i]) for i in range(self.nrsig)])
          self.diagerr = diagerr
        except:
          diagerr = np.zeros(self.nrsig) + 1.
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

#  def __sub__(self,other):
#    if(self.checkbinning(other) == 0):
#      wpnew = copy.deepcopy(self)
#      wpnew.wp = other.wp - self.wp
#    else:
#      print 'add fail!  binning misaligned'
#      wpnew = None
#    return wpnew

  def printwp(self,outfname):
    ofp = open(outfname,'w')
    if(self.DorT == 1): ##this means has no covariance matrix.
      for i in range(len(self.wp)):
        ofp.write('%e %e\n' % (self.rsig[i],self.wp[i]))
    else:
      for i in range(len(self.wp)):
        ofp.write('%e %e %e\n' % (self.rsig[i],self.wp[i],self.diagerr[i]))
    ofp.close()

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
    elif self.DorT == 0 and other.DorT == 0:
      return ((diff*self.icov) * (diff.T))[0,0]
    else:
      return 0. ## can't compute chi2 without an inverse covariance matrix.

  def chi2hack(self,other,otherstart,otherend):
    """
    compute chi2 of self with other.  self is the one with the inverse cov.
    """
    diff = np.matrix(self.wp-other.wp[otherstart:otherend+1])
    if self.DorT == 0 and other.DorT == 1:
      return ((diff*self.icov) * (diff.T))[0,0]
    elif self.DorT == 1 and other.DorT == 0:
      return ((diff*other.icov) * (diff.T))[0,0]
    else:
      return 0. ## can't compute chi2 without an inverse covariance matrix.

  def chi2info(self,other):
    """
    For non-diagonal cov, this helps show you which differences between model and theory contribute the most
    to chi2.
    """
    if isinstance(other,wp):
      diff = np.matrix(self.wp-other.wp)
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
        diff = np.matrix(self.wp - other)
        myicov = self.icov
      except:
        return 0.
    mysum = 0.
    for i in range(len(self.wp)):
      print i,diff[0,i]*((diff*myicov)[0,i])
      mysum += diff[0,i]*((diff*myicov)[0,i])
    assert np.fabs(mysum - self.chi2(other)) < 0.001
    return mysum


def wpfromDDoDR(fbase,dfacr=1,periodicopt=0,DRfacinfo=None,rpimax=80.,testing=0,icovfname=None):
  """
  This function computes wp from DRopt[1-2] files.
  DRfacinfo is used to pass [DRfac, fixRR] (dont get it from file headers in this case)
  That feature is necessary for computing covariance matrices from bootstrap, where that factor
  should be fixed by total N and/or S, it does not vary from bootstrap region to region.
  """

  if(dfacr != 1):
    print 'dfacr not 1 is not coded up yet!'
    sys.exit(1)

  if(periodicopt == 1): # periodic box.
    print 'did not code up periodic yet!'
    sys.exit(1)

  else:
    fDD = fbase+'.DRopt1'
    fDR = fbase+'.DRopt2'
    rpg, rpig, DDg = np.loadtxt(fDD,skiprows=3,unpack=True)
    rpg, rpig, DRg = np.loadtxt(fDR,skiprows=3,unpack=True)

    if DRfacinfo is None:
      DRfac = ximisc.getDRfactoronly(fbase)
    else:
      try:
        DRfac = DRfacinfo[0]
      except:
        DRfac = DRfacinfo

    print 'using DRfac',DRfac

    nrpibins = len(np.where(rpg == rpg[0])[0])
    nrpbins = len(np.where(rpig == rpig[0])[0])
    if nrpibins*nrpbins != len(DDg):
      return None

    xi = (2.*DDg/DRg/DRfac - 1.)

    rpi1d = rpig.reshape(nrpbins,nrpibins)[0]
    ## check that it's linearly spaced, get drpi
    drpi = rpi1d[1]-rpi1d[0]
    chk = rpi1d[1:] - rpi1d[:-1] - drpi
    aaa = np.where(np.fabs(chk) > 1.0e-5)[0]
    if(len(aaa) > 0):
      print 'error!'
      sys.exit(1)
    assert np.fabs(rpi1d[0]-drpi*0.5) < 1e-3
    rp1d = rpg.reshape(nrpbins,nrpibins)[:,0]

    mywp = np.zeros(len(rp1d),dtype='float')
    wpi = 0
    for rpval in rp1d:
      if(testing==1):
        print 'inside testing'
        xx = np.where(rpg == rpval)[0]
        assert len(xx) == nrpibins
        xicurr = xi[xx]
        assert (rpig[xx] == rpi1d).all()
        assert np.fabs(rpi1d[0]-drpi*0.5) < 1e-3

      ## ok, we're sure rpigrid is sane.
      xx = np.where((rpg == rpval) & (rpig < rpimax))[0]
      mywp[wpi] = (xi[xx]).sum()*drpi*2.
      wpi += 1

    return wp(wpfname=fDR,icovfname=icovfname,rpwplist=[rp1d,mywp])
 

def wpfromDR(fbase,dfacr=1,periodicopt=0,DRfacinfo=None,rpimax=80.,testing=0,icovfname=None,wpstart=-1,wpend=-1,xiopt=0):
  """
  This function computes wp from DRopt[1-3] files.
  DRfacinfo is used to pass [DRfac, fixRR] (dont get it from file headers in this case)
  That feature is necessary for computing covariance matrices from bootstrap, where that factor
  should be fixed by total N and/or S, it does not vary from bootstrap region to region.
  xiopt added to input Hong's xi values rather than DR counts.
  """
  if(dfacr != 1):
    print 'dfacr not 1 is not coded up yet!'
    sys.exit(1)

  if(periodicopt == 1): # periodic box.
    print 'did not code up periodic yet!'
    sys.exit(1)
  elif xiopt == 1:
    rpg, rpig, xi, xigerr = np.loadtxt(fbase,unpack=True)
    nrpibins = len(np.where(rpg == rpg[0])[0])
    nrpbins = len(np.where(rpig == rpig[0])[0])
    if nrpibins*nrpbins != len(xi):
      return None
    fDR = fbase

  else:
    fDD = fbase+'.DRopt1'
    fDR = fbase+'.DRopt2'
    fRR = fbase+'.DRopt3'
    rpg, rpig, DDg = np.loadtxt(fDD,skiprows=3,unpack=True)
    rpg, rpig, DRg = np.loadtxt(fDR,skiprows=3,unpack=True)
    rpg, rpig, RRg = np.loadtxt(fRR,skiprows=3,unpack=True)

    if DRfacinfo is None:
      DRfac, fixRR = ximisc.getDRfactors(fbase)
    else:
      DRfac = DRfacinfo[0]
      fixRR = DRfacinfo[1]

    nrpibins = len(np.where(rpg == rpg[0])[0])
    nrpbins = len(np.where(rpig == rpig[0])[0])
    if nrpibins*nrpbins != len(DDg):
      return None

    xi = ((DDg-DRg*DRfac)/RRg/DRfac**2 + 1.)

  rpi1d = rpig.reshape(nrpbins,nrpibins)[0]
  ## check that it's linearly spaced, get drpi
  drpi = rpi1d[1]-rpi1d[0]
  chk = rpi1d[1:] - rpi1d[:-1] - drpi
  aaa = np.where(np.fabs(chk) > 1.0e-5)[0]
  if(len(aaa) > 0):
    print 'error!'
    sys.exit(1)
  assert np.fabs(rpi1d[0]-drpi*0.5) < 1e-3
  rp1d = rpg.reshape(nrpbins,nrpibins)[:,0]

  mywp = np.zeros(len(rp1d),dtype='float')
  wpi = 0

  for rpval in rp1d:
    if(testing==1):
      print 'inside testing'
      xx = np.where(rpg == rpval)[0]
      assert len(xx) == nrpibins
      xicurr = xi[xx]
      assert (rpig[xx] == rpi1d).all()
      assert np.fabs(rpi1d[0]-drpi*0.5) < 1e-3

    ## ok, we're sure rpigrid is sane.
    xx = np.where((rpg == rpval) & (rpig < rpimax))[0]
    mywp[wpi] = (xi[xx]).sum()*drpi*2.
    wpi += 1

  return wp(wpfname=fDR,icovfname=icovfname,rpwplist=[rp1d,mywp],wpstart=wpstart,wpend=wpend)



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


