import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ximisc
import xiell

class xismu:
  def __init__(self,smufname=None,icovfname=None,Npopt=0,smuxilist=[],smuvNplist=[]):
    """
    Read in a correlation function measurement of the form
    s, mu, xi.
    If Npopt = 1, assumes you have an "Np" file, which has the format
    s, mu, vol, Npairs.  xi will be computed from pair counts.
    There are methods to rebin and plot.
    """

    ## make these exist but none so I can check for their presence
    ## when I want to use them.

    self.Np = None
    self.vol = None

#    try:
    if(0==0):
      if(Npopt == 1 and smuvNplist != []):
        self.s = smuvNplist[0]
        self.mu = smuvNplist[1]
        self.vol = smuvNplist[2]
        self.Np = smuvNplist[3]
        self.xi = self.Np/self.vol-1.
      elif(Npopt == 1 and smufname is not None):
        self.Npsmufname = [smufname]
        self.s, self.mu, self.vol, self.Np = np.loadtxt(smufname,unpack=True)
        self.xi = self.Np/self.vol-1.
      elif(Npopt == 0 and smufname is not None):
        self.xismufname = [smufname]
        self.s, self.mu, self.xi = np.loadtxt(smufname,unpack=True)
      elif smuxilist != []:
        self.s = smuxilist[0]
        self.mu = smuxilist[1]
        self.xi = smuxilist[2]
      else:  ## need one of those input cases to be true.
        assert 0==1

      self.ntot = len(self.s)
      ## figure out binning yourself.
      self.nmu = len(np.where(self.s == self.s[0])[0])
      self.ns = self.ntot/self.nmu
      ## no remainder, I got the binning right.
      assert (self.ntot % self.nmu) == 0
      self.s1d = self.s.reshape(self.ns, self.nmu)[:,0]
      self.mu1d = self.mu.reshape(self.ns, self.nmu)[0,:]

      ## compute linear spacings.
      dmu = (self.mu1d[1:]-self.mu1d[:-1]).mean()
      ds = (self.s1d[1:]-self.s1d[:-1]).mean()
      if (np.fabs(self.mu1d[1:]-self.mu1d[:-1] - dmu) > 0.0001*dmu).any():
        ## not linear binning!
        ## this shouldn't be the case for a s,mu grid.
        self.dmu = 0.
        assert self.dmu > 0.
      else:
        self.dmu = dmu 

      if (np.fabs(self.s1d[1:]-self.s1d[:-1] - ds) > 0.0001*ds).any():
        ## not linear binning!
        self.ds = 0.
      else:
        self.ds = ds
        self.logsopt = 0

      ## check for log binning.
      dlogs = (np.log(self.s1d[1:]/self.s1d[:-1])).mean()
      if (np.fabs(np.log(self.s1d[1:])-np.log(self.s1d[:-1]) - dlogs) < 0.0001*dlogs).all():
        self.logsopt = 1
        self.dlogs = dlogs

#    except:
    else:
      print 'bad inputs',smufname, smuxilist
      self = None

  def simplexiell(self):
    """
    Computes multipoles in the simplest way -- 
    Direct summation over the bins times Legendre polynomials.
    No interpolation or correction for the fact that xi_{2,4,...} may not
    be 0 for xi=0.
    """
    myxiell = np.zeros([3,self.ns],dtype=float)
    mysvec = np.zeros([3,self.ns],dtype=float)
    for i in range(self.ns):
      i1 = i*self.nmu
      i2 = (i+1)*self.nmu
      mymu = self.mu[i1:i2]
      ## sanity check that you picked out a correct and full range of mu.
      assert (mymu[1:] - mymu[:-1] > 0.).all()
      dmu = (mymu[1:] - mymu[:-1]).mean()
      assert np.fabs(mymu[0]-dmu*0.5) < 2.0e-4
      assert np.fabs(mymu[-1]+dmu*0.5-1.) < 2.0e-4
      for ell in [0,2,4]:
        myxiell[ell/2,i] = (self.xi[i1:i2]*ximisc.legendre(ell,mymu)).sum()*(2.*ell+1.)/float(self.nmu)
        mysvec[ell/2,i] = self.s1d[i]
        ## tmp!
        #print '%e %e\n' % (mysvec[ell/2,i], myxiell[ell/2,i])
      ## create a xiell object.
    return xiell.xiell(ellmax=4,sxilist=[mysvec,myxiell])



  def printsmuxi(self,outfname):
    """ Print a file with columns s, mu, xi."""
    ofp = open(outfname,'w')
    for i in range(self.ntot):
      ofp.write('%e %e %e\n' % (self.s[i], self.mu[i], self.xi[i]))
    ofp.close()

  def rebin(self,dfacs,dfacmu):
    """
    Rebin by factor dfacs/dfacmu in s, mu.
    Creates a new xismu object and returns it.
    """

    if self is None:
      return None
    ## can't do this calculation without pair counts/volumes.
    if(self.Np is None or self.vol is None):
      return None

    ## make sure the binning works out ok; its ok if it's not even for ns
    if self.nmu % dfacmu != 0:
      return None 

    nmudown = self.nmu/dfacmu
    ristart = 0
    nsdown = (self.ns)/dfacs
    #### don't need these.
    sdown = ximisc.downsample1d(self.s1d,dfacs)
    assert len(sdown) == nsdown
    mudown = ximisc.downsample1d(self.mu1d,dfacmu)
    ## this is more general than equal bins.
    xinew = np.zeros([nsdown,nmudown])


    for i in range(nsdown):
      for ishort in range(dfacs):
        i1 = (ristart+ishort)*self.nmu
        i2 = (ristart+ishort+1)*self.nmu
        if(ishort == 0):
          myvol = ximisc.downsample1dsum(self.vol[i1:i2],dfacmu)
          myNp = ximisc.downsample1dsum(self.Np[i1:i2],dfacmu)
        else:
          myvol = myvol + ximisc.downsample1dsum(self.vol[i1:i2],dfacmu)
          myNp = myNp + ximisc.downsample1dsum(self.Np[i1:i2],dfacmu)
  
      xinew[i,:] = myNp/myvol-1.
      ristart += dfacs

    s2d, mu2d = ximisc.one2twod(sdown,mudown)
    #x = copy.deepcopy(self)
    #xismu.__init__(x,smuxilist=[s2d.flatten(),mu2d.flatten(),xinew])
    #return x
    return xismu(smuxilist=[s2d.flatten(),mu2d.flatten(),xinew.flatten()])

  def ximubin(self,mubin):
    return(self.xi.reshape(self.ns,self.nmu)[:,mubin])

  def addcurve(self,ax,mubin,color='k',spow=2,fmt=None,lbl=None):
    """
    adds a curve on a xis(mubin) plot. specify which mubin to plot.
    """
    mys = self.s1d
    plt.plot(mys,self.ximubin(mubin)*mys**spow,color=color,label=lbl)

  def makexismuplot(self,panelopt=0,color='k',span=None,logxopt=1,logyopt=0,spow=2,fmt=None,lbl=None):
    """
    Color should be a single color or a list of length self.nmu
    """


    try:
      cval = color[self.nmu-1]
    except:
      tmp = []
      for i in range(self.nmu):
        tmp.append(color)
      color = tmp

    if lbl is not None:
      try:
        lval = lbl[self.nmu-1]
      except:
        tmp = []
        for i in range(self.nmu):
          tmp.append(lbl)
        lbl = tmp


    if(self.nmu > 5 and panelopt == 1):
      print 'too many panels to make, switching to panelopt = 0.'
      panelopt = 0

    if(panelopt == 1):
      if self.nmu == 5:
        nx = 3
        ny = 2
      if self.nmu == 4:
        nx = 2
        ny = 2
      if self.nmu == 3:
        nx = 1
        ny = 3
      if self.nmu == 2:
        nx = 1
        ny = 2
   
      self.fig = plt.figure(figsize=[5*nx,5*ny])
      axlist = []
      for i in range(self.nmu):
        axlist.append(self.fig.add_subplot(nx,ny,i+1))
        if lbl is not None:
          self.addcurve(axlist[i],i,color=color[i],spow=spow,lbl=lbl[i],fmt=fmt)
        else:
          self.addcurve(axlist[i],i,color=color[i],spow=spow,lbl=lbl,fmt=fmt)

    else: ## put everything on one
      self.fig = plt.figure(figsize=[6,6])
      axlist = []
      axlist.append(self.fig.add_subplot(1,1,1))
      for i in range(self.nmu):
        if lbl is not None:
          self.addcurve(axlist[0],i,color=color[i],spow=spow,lbl=lbl[i],fmt=fmt)
        else:
          self.addcurve(axlist[0],i,color=color[i],spow=spow,lbl=lbl,fmt=fmt)

    ###### nevermind, let's just let it be auto-set here if nothing input.
    #if span is None:
    #  mys = self.s1d
    #  span = [mys.min()*0.9,mys.max()*1.1,((mys)**spow*(self.ximubin(mubin))).min()*0.9,((mys)**spow*(self.ximubin(mubin))).max()*1.1]

    for aa in axlist:
      if(logxopt == 1):
        aa.set_xscale('log')
      if(logyopt == 1):
        aa.set_yscale('log')

      if span is not None:
        aa.axis(span)
      aa.set_xlabel(r'$s \, [h^{-1} {\rm Mpc}]$',fontsize=20)
      if(np.fabs(spow) > 0.01):
        aa.set_ylabel(r'$s^{%.1f} \xi_{\mu}(s)$' % (spow),fontsize=20)
      else:
        aa.set_ylabel(r'$\xi_{\mu}(s)$',fontsize=20)

    return axlist



def initavg(ffmt,zspacelist,simlist):
  cnt = 0
  for s in simlist:
    for z in zspacelist:
      ff = ffmt % (s,z)
      aa = np.loadtxt(ff)
      if(cnt==0):
        aatot = copy.deepcopy(aa)
        svec = aa[:,0]
        muvec = aa[:,1]
      else:
        aatot = aatot + aa
## make sure all the files have the same binning.
        assert (np.fabs(svec - aa[:,0]) < 2.0e-4).all()
        assert (np.fabs(muvec - aa[:,1]) < 2.0e-4).all()
      cnt += 1
  aa = aatot/float(cnt)
  (nr, nc) = aa.shape
  if nc == 4:
    myl = [aa[:,0],aa[:,1],aa[:,2],aa[:,3]]
    return xismu(Npopt=1,smuvNplist=myl)
  if nc == 3:
    myl = [aa[:,0],aa[:,1],aa[:,2]]
    return xismu(Npopt=0,smuxilist=myl)

  ## don't know what to do with a different number of columns.
  return None


if __name__ == '__main__':
  ##tmp testing.
  aa=xismu(smufname="/Users/bareid/work/montserratdata/forlilev3/A03_0.6452.halos.zspace0.Np.nbins9000.bin1sims_Mmin12.182Mmax12.483",Npopt=1)
  arebin=aa.rebin(1,1)
  #arebin.printsmuxi("tmpodebug.dat")
  myxiell=arebin.simplexiell()
  myxiell.printxiellshort("xiellchk.dat")

  ## works!!
#bash-3.2$ diff xiellchk.dat /Users/bareid/work/montserratdata/forlilev3/A03_0.6452.halos.zspace0.xiell.nbins9000.bin1sims_Mmin12.182Mmax12.483
#bash-3.2$ 

  ## test the averaging thing.
  ffmt="/Users/bareid/work/montserratdata/forlilev3/A%02d_0.6452.halos.zspace%d.Np.nbins9000.bin1sims_Mmin12.182Mmax12.483"
  xx=initavg(ffmt,[0,1,2],[0,2,3,4,5])
  xx.printsmuxi("avgmbin1.smuxi")
  #def printsmuxi(self,outfname):

  xrebin=xx.rebin(1,50)
  xrebin.printsmuxi("avgmbin1down.smuxi")
  xrebin.makexismuplot(panelopt=0,logxopt=0)
  xrebin.fig.savefig("testoyo.pdf")
  #def makexismuplot(self,panelopt=0,color='k',span=None,logxopt=1,logyopt=0,spow=1,fmt=None,lbl=None):
  
  
