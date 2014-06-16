import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
#import ximisc
#import cov
import re
import scipy.special
import BRutils
import scipy.stats
import time
  
def normGauss(v,sig2):
  return np.exp(-0.5*v**2/sig2)/np.sqrt(2.*np.pi*sig2)

class massfxn:
  def __init__(self,m=None,Nofm=None,mffname=None,Lbox=1.):
    """
    To rescale counts to number densities per ln M (not ln10!), fill an a value for Lbox.  Leave as one if you're inputing number densities already and nothing will be rescaled.  Else scale by Lbox**3 * [dm or dlgm = dlg10m * np.log(10.)]
    """
    if m is not None:
      self.m = m
      self.Nofm
      self.log10mcen = np.log10(self.m) ## this mass defn matches the class halocat definition.
      self.logopt = -1

    else:
      if mffname is not None:
        myfmt = 0
        ## current format from precompute.
        ifp = open(mffname,'r')
        line = ifp.readline()
        if(re.search('# dlg10m:',line)):
          self.dlg10m = float(line.split(':')[1].strip())
          myfmt += 1
        line = ifp.readline()
        if(re.search('# minlg10m:',line)):
          self.lg10mmin = float(line.split(':')[1].strip())
          myfmt += 1
        ifp.close()
        if myfmt == 2: ## it's lg10m masses.
          self.lg10mcen, self.Nofm = np.loadtxt(mffname,usecols=[0,1],unpack=True)
          self.m = 10**self.lg10mcen
          assert ((self.lg10mcen - self.lg10mmin - 0.5*self.dlg10m)/self.dlg10m - np.arange(0,len(self.lg10mcen),1) < 0.001).all()
        else:
          self.m, self.Nofm = np.loadtxt(mffname,usecols=[0,1],unpack=True)
          self.log10mcen = np.log10(self.m)


      else:
        self = None
        return None

    dm = (self.m[1:] - self.m[:-1]).mean()
    dlg10m = (np.log10(self.m[1:]/self.m[:-1])).mean()
    if(np.fabs((self.m[1:] - self.m[:-1])/dm -1.)  < 0.001).all():
      self.logopt = 0
      self.dm = dm
      if(np.fabs(Lbox - 1.) > 2.0e-6):
        self.Nofm = self.Nofm/Lbox**3/self.dm
    if(np.fabs(np.log10(self.m[1:]/self.m[:-1])/dlg10m - 1.) < 0.001).all():
      self.logopt = 1
      self.dlg10m = dlg10m
      if(np.fabs(Lbox - 1.) > 2.0e-6):
        self.Nofm = self.Nofm/Lbox**3/(self.dlg10m*np.log(10.))

  def addmassfxncurve(self,ax=None,color='k',lbl=None):
    """copied from the elt of sims.py"""
    if ax is None:
      plt.plot(10**(self.lg10mcen),self.Nofm,color=color,label=lbl)
    else:
      ax.plot(10**(self.lg10mcen),self.Nofm,color=color,label=lbl)

  def massfxnplot(self,color='k',lbl=None):
    """copied from the elt of sims.py"""
    f=plt.figure()
    ax = f.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    self.addmassfxncurve(ax,color,lbl)
    return f, ax




class hod():
  def fitvdisppowerlaw(self, Mmin=4.e12, Mmax=5.e14):
#    print 'hi beth',len(self.vdisp), len(self.mf.m)
    if 0==0:
      xx = np.where((self.vdisp > 0) & (self.mf.m >= Mmin) & (self.mf.m <= Mmax))[0]
#      print len(xx)
      slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log(self.mf.m[xx]/1.0e13),np.log(self.vdisp[xx]))
      self.vdispnorm = np.exp(intercept)
      self.vdisppow = slope
    else:
      print 'could not fit vdisp!'
      

  ## default HOD is current beth best fit to xiell's.
  def __init__(self,Mmin=1.1460142e+13,M1=1.57663077e+14,alpha=1.298719e+00,Mcut=3.238620e+11,sigmalogM=3.40466e-01,m=None,Nofm=None,mffname=None,Lbox=677.7,cenopt=2,satopt=2,vdispfname=None):

    if(mffname is not None):
      self.mf = massfxn(mffname=mffname,Lbox=Lbox,m=m,Nofm=Nofm)

    if mffname is not None and vdispfname is not None:
      m1 = np.loadtxt(mffname,unpack=True,usecols=[0])
      m2 = np.loadtxt(vdispfname,unpack=True,usecols=[0])
      assert (m1 == m2).all()


    if(vdispfname is not None):
      if mffname is None:
        self.mf = massfxn(mffname=vdispfname,Lbox=Lbox,m=m,Nofm=Nofm)
      self.vdisp = np.loadtxt(vdispfname,unpack=True,usecols=[4])
      ## fit a smooth version.
      self.fitvdisppowerlaw()
      
 
    self.Mmin = Mmin
    self.M1 = M1
    self.alpha = alpha
    self.Mcut = Mcut
    self.sigmalogM = sigmalogM
    self.cenopt = cenopt
    self.satopt = satopt


  def fsat(self):
    try:
      fsat = (self.mf.Nofm*self.Ncen(self.mf.m)*(self.Nsat(self.mf.m))).sum()/(self.mf.Nofm*self.Ncen(self.mf.m)*(1+self.Nsat(self.mf.m))).sum()
    except:
      fsat = -1.
    return fsat
  
  def nbar(self):
    try:
      nbar = ((self.mf.dlg10m*np.log(10.)) * self.mf.Nofm*self.Ncen(self.mf.m)*(1+self.Nsat(self.mf.m))).sum()
    except:
      nbar = -1.
    return nbar

  def vdispsmooth(self,m):
    return self.vdispnorm*(m/1.0e13)**self.vdisppow    

  def vrmssmooth(self,m):
    return (self.vdispsmooth(m))**0.5

  def Ncen(self,M):
    """
    works for scalar or vector inputs M.
    """
    try:
      len(M)
      Mvec = M
      scalarflag = 0
    except:
      Mvec = np.array([M])
      scalarflag = 1
    result = np.zeros(len(Mvec))

    if(self.cenopt == 0):
      result[:] = 1.0


    if(self.cenopt == 1):
      xx = np.where(Mvec < self.Mmin)[0]
      yy = np.where(Mvec >= self.Mmin)[0]
      result[xx] = 0.0
      result[yy] = 1.0

    if(self.cenopt == 2 or self.cenopt == 5):
      if(self.sigmalogM < 1.e-4):
        xx = np.where(Mvec < self.Mmin)[0]
        yy = np.where(Mvec >= self.Mmin)[0]
        result[xx] = 0.0
        result[yy] = 1.0
      else:
        result = 0.5*(1.+scipy.special.erf(np.log10(Mvec/self.Mmin)/self.sigmalogM))

    if(self.cenopt == 15):
      xx = np.where((Mvec >= self.Mmin) & (Mvec < self.M1))[0]
      result[xx] = 1.
    if scalarflag == 1:
      return result[0]
    else:
      return result

  def Nsat(self,M):
    """
    works for scalar or vector inputs M.
    """
    try:
      len(M)
      Mvec = M
      scalarflag = 0
    except:
      Mvec = np.array([M])
      scalarflag = 1
    result = np.zeros(len(Mvec))

    if self.satopt == 0:
      ##
      pass ## already result is 0's.
    if self.satopt == 1:
      xx = np.where(Mvec < self.Mmin)[0]
      yy = np.where(Mvec >= self.Mmin)[0]
      result[xx] = 0.
      result[yy] = pow(Mvec[yy]/self.M1,self.alpha)
    if self.satopt == 2:
      xx = np.where(Mvec < self.Mcut)[0]
      yy = np.where(Mvec >= self.Mcut)[0]
      result[xx] = 0.
      result[yy] = pow((Mvec[yy]-self.Mcut)/self.M1,self.alpha)

    if scalarflag == 1:
      return result[0]
    else:
      return result

  def plothod(self,ax=None, Mlow=1.e11,Mhigh=5e15,dlogM=0.01,popt=2,color='k',MFwgt=False,rebinsize=1,renormalize=False):
    """
    popt = 0: cen only
    popt = 1: sat only
    popt = 2: tot
    """

    if MFwgt == True:
      Mlist = self.mf.m
      wgt = self.mf.Nofm
      if rebinsize > 1:
        Mlist = np.exp(BRutils.rebin(np.log(self.mf.m),rebinsize=rebinsize))
        tmpo = np.log10(Mlist[1:]/Mlist[:-1]).mean()
        assert (np.fabs(np.log10(Mlist[1:]) - np.log10(Mlist[:-1]) - tmpo) < 0.001*tmpo).all()
        Nrebin = BRutils.rebin(self.mf.Nofm,rebinsize=rebinsize)
        wgt = Nrebin
        ## need to rescale by dlog10M fac??
    else:
      Mlist = np.exp(np.arange(np.log(Mlow),np.log(Mhigh)+0.5*dlogM,dlogM))
      wgt = np.zeros(len(Mlist))+1.

    if ax is None:
      ff = plt.figure(figsize=[6,6])
      ax = ff.add_subplot(1,1,1)
    else:
      ff = None
    if popt == 0:
      myy = wgt*self.Ncen(Mlist)
    if popt == 1:
      myy = wgt*self.Ncen(Mlist)*(0.+self.Nsat(Mlist))
    if popt == 2:
      myy = wgt*self.Ncen(Mlist)*(1.+self.Nsat(Mlist))
    if renormalize == True:
      myy = myy/float(myy.sum())


    ax.plot(Mlist,myy,color)
    ax.set_xscale('log')
    return ff, ax


  def fillpofv(self,dv=0.05,cenvfrac=0.3):
    vmax = self.vrmssmooth(self.mf.m.max())*3.
    self.vlist = np.arange(0.,vmax,dv)
    self.pvlist = np.zeros(len(self.vlist))
    self.cenvfrac=cenvfrac
    print 'need to do',len(self.vlist),'times',len(self.mf.m)
    ntot = 0.
    for mval, Nofmval in zip(self.mf.m,self.mf.Nofm):
      ntot += Nofmval*self.Ncen(mval)*(1.+self.Nsat(mval))

    for ii in range(len(self.vlist)):
      if ii % 100 == 0:
        print ii
      vv = self.vlist[ii]
      ## sum over mass fxn.
      for mval, Nofmval in zip(self.mf.m,self.mf.Nofm):
#        print dv, Nofmval, self.Ncen(mval), self.Nsat(mval), normGauss(vv,cenvfrac**2*self.vdispsmooth(mval)), dv*Nofmval*self.Ncen(mval)*normGauss(vv,cenvfrac**2*self.vdispsmooth(mval)), dv*Nofmval*self.Ncen(mval)*self.Nsat(mval)*normGauss(vv,self.vdispsmooth(mval))
        self.pvlist[ii] += dv*Nofmval*self.Ncen(mval)*normGauss(vv,cenvfrac**2*self.vdispsmooth(mval)) + \
                    dv*Nofmval*self.Ncen(mval)*self.Nsat(mval)*normGauss(vv,self.vdispsmooth(mval))
      self.pvlist[ii] = self.pvlist[ii]/ntot*2. # 2. for negative values.
#      print 'ack',self.pvlist[ii] 

  def sig2hod(self,cenvfrac=0.3,ihvscale=1.0):
    """
    Returns total sig2 given a cenvfrac.
    """
    return (self.mf.Nofm*self.Ncen(self.mf.m) * \
           (cenvfrac**2 + ihvscale**2 * self.Nsat(self.mf.m))*self.vdispsmooth(self.mf.m)).sum()/ \
           (self.mf.Nofm*self.Ncen(self.mf.m)*(1.+self.Nsat(self.mf.m))).sum()

  def sig2hodcum(self,cenvfrac=0.3,ihvscale=1.0):
    """
    Returns total sig2 given a cenvfrac.
    """
    return (self.mf.Nofm*self.Ncen(self.mf.m) * \
           (cenvfrac**2 + ihvscale**2 * self.Nsat(self.mf.m))*self.vdispsmooth(self.mf.m)).cumsum()/ \
           (self.mf.Nofm*self.Ncen(self.mf.m)*(1.+self.Nsat(self.mf.m))).sum()

def hodfromchain(celt,whichbox=1,cenopt=2,satopt=2,bosswdir='/home/howdiedoo/boss/',whichPB=0,COMVopt=0):
  """
  similar structure to hodutils.py/runchainmodel
  whichbox = 0 for L0, 1 for MWhires, 2 for PB
  """
  logopt = 0
  if celt['M_min'] < 1e10:
    logopt = 1

  if logopt == 1:
    Mmin = 10**(celt['M_min'])
    Mcut = 10**(celt['M_cut'])
    M1 = 10**(celt['M1'])
  else:
    Mmin = (celt['M_min'])
    Mcut = (celt['M_cut'])
    M1 = (celt['M1'])

  if whichbox == 0:
    Lbox = 2750.0
    mffname = bosswdir + '/bethalexieprecpinterp/precompute_files/precomputenewvelcombo_HV1.000_IHV1.000_nmh_nsat.txt'
    vdispfname = None
  if whichbox == 1:
    Lbox = 677.7
    ## doesn't exist on lappy
    #mffname = bosswdir + '/bethalexieprecpinterpNOFOF/precompute_files/precompute_MWhiresplanck_nmh_nsat.txt'
    mffname = bosswdir.split('/boss/')[0]+'/SOmaster/MWhiresplanck.halomembers.vdisp'
    vdispfname = bosswdir.split('/boss/')[0]+'/SOmaster/MWhiresplanck.halomembers.vdisp'
    if COMVopt == 1:
      vdispfname = bosswdir.split('/boss/')[0]+'/SOmaster/MWhiresplanckCOMV.halomembers.vdisp'
  if whichbox == 2:
    Lbox = 1380.0
    mffname = bosswdir + 'bethalexieprecpinterpNOFOFxiellwpwcov/precompute_files/precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000_nmh_nsat.txt' % whichPB
    vdispfname = bosswdir.split('/boss/')[0]+'/SOmaster/PB%02d_0.6452.halomembers.vdisp' % whichPB
    if COMVopt == 1:
      vdispfname = bosswdir.split('/boss/')[0]+'/SOmaster/PB%02dCOMV_0.6452.halomembers.vdisp' % whichPB

  hh = hod(Mmin = Mmin, M1 = M1, Mcut = Mcut, sigmalogM= celt['sigma_logM'], alpha = celt['alpha'],mffname=mffname,Lbox=Lbox,cenopt=cenopt,satopt=satopt,vdispfname=vdispfname)
  return hh

def comparevcats(whichbox=1,bosswdir='/home/howdiedoo/boss/',whichPB=0,sanitycheck=1,msplits=None):
  if whichbox == 1:
    fCOMV = bosswdir.split('boss/')[0] + 'SOmaster/MWhiresplanckCOMV.halos'
    fnew = bosswdir.split('boss/')[0] + 'SOmaster/MWhiresplanck.halos'
    fvCOMV = bosswdir.split('boss/')[0] + 'SOmaster/MWhiresplanckCOMV.halomembers.vdisp'
    fvnew = bosswdir.split('boss/')[0] + 'SOmaster/MWhiresplanck.halomembers.vdisp'
    Lbox = 677.7

  if whichbox == 2:
    Lbox = 1380.0
    assert whichPB == 0 ## only one we ran COMV for.
    fCOMV = bosswdir.split('/boss/')[0]+'/SOmaster/PB%02dCOMV_0.6452.halos' % whichPB
    fnew = bosswdir.split('/boss/')[0]+'/SOmaster/PB%02d_0.6452.halos' % whichPB
    fvCOMV = bosswdir.split('/boss/')[0]+'/SOmaster/PB%02dCOMV_0.6452.halomembers.vdisp' % whichPB
    fvnew = bosswdir.split('/boss/')[0]+'/SOmaster/PB%02d_0.6452.halomembers.vdisp' % whichPB

  #ifp1 = open(fCOMV,'r')
  #ifp2 = open(fnew,'r')

  if sanitycheck == 1:
    a1 = np.loadtxt(fCOMV,unpack=True,usecols=[0,1,2,7])
    a2 = np.loadtxt(fnew,unpack=True,usecols=[0,1,2,7])
    assert(a1 == a2).all()
    print 'passed assert!'
    del a1
    del a2

  ## positions and masses are the same, just read in velocities and masses.
  t1=time.time()
  m = np.loadtxt(fCOMV,unpack=True,usecols=[7])
  vCOMV = np.loadtxt(fCOMV,unpack=True,usecols=[3,4,5])
  vnew = np.loadtxt(fnew,unpack=True,usecols=[3,4,5])
  xx = m.argsort()
  m = m[xx]
  for ii in [0,1,2]:
    vCOMV[ii,:] = vCOMV[ii,xx]
    vnew[ii,:] = vnew[ii,xx]
  t2=time.time()
  print 'read and sort took:',t2-t1

  if msplits is None:
    ## doesn't work for PB
    #msplits = np.array([10**m.min()*0.99,1e12,5e12,8e12,1e13,2e13,5e13,1e14,5e14,10**m.max()*1.01])
    msplits = np.array([1.27e12,5e12,8e12,1e13,2e13,5e13,1e14,5e14,10**m.max()*1.01])
    msplits = np.log10(msplits)
  illist=[]
  ihlist=[]
  mmedlist=[]
  for ml, mh in zip(msplits[:-1], msplits[1:]):
    tmp = np.where((m >= ml) & (m < mh))[0]
    illist.append(tmp.min())
    ihlist.append(tmp.max()+1)
    mmedlist.append((m[tmp])[len(tmp)/2])

  nmchunk = len(mmedlist)
  nstats = 7
  s = np.zeros([nmchunk,nstats])
  for mi in range(nmchunk):
    #istart = nhalos/nmchunk*mi
    #iend = nhalos/nmchunk*(mi + 1)+1
    #if mi == nmchunk - 1: iend = nhalos
    istart = illist[mi]
    iend = ihlist[mi]
    s[mi,0] = mmedlist[mi] #m[(istart+iend)/2] ## close enough.
    for ii in [0,1,2]:
      s[mi,1] += ((vCOMV[ii,istart:iend] - vnew[ii,istart:iend])**2).sum()
      s[mi,2] += ((vCOMV[ii,istart:iend])**2).sum()
      s[mi,3] += ((vnew[ii,istart:iend])**2).sum()
      s[mi,4] += ((vnew[ii,istart:iend])*(vCOMV[ii,istart:iend])).sum()
      s[mi,5] += ((vCOMV[ii,istart:iend] - vnew[ii,istart:iend])*(vCOMV[ii,istart:iend])).sum()
      s[mi,6] += ((vCOMV[ii,istart:iend] - vnew[ii,istart:iend])*(vnew[ii,istart:iend])).sum()
    for jj in range(1,7):
      s[mi,jj] = s[mi,jj]/(3.*(iend-istart))
      if jj <= 3:
        s[mi,jj] = s[mi,jj]**0.5*Lbox
    s[mi,4] = s[mi,4]*Lbox**2/(s[mi,2]*s[mi,3])
    s[mi,5] = s[mi,5]*Lbox**2/(s[mi,1]*s[mi,2])
    s[mi,6] = s[mi,6]*Lbox**2/(s[mi,1]*s[mi,3])

  ## return set of interesting statistics about velocity differences.
  return s
  


  




if __name__ == '__main__':
  pass


