import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import scipy.optimize
import os
import re

from scipy.integrate import quad
import scipy.interpolate as interp

## stuff needed to describe neutrino energy density evolution.
global finterp
global Ibigy
global yi, fyi
global onuh2const
global setup 
global dm21sqr, dm31sqr
global Smnumin  #array of size 2. indexed by NorI

global rhocrit

setup = 0
onuh2const = 7./8.*(4./11.)**(4./3.)
print 'onuhu2 conversion!',onuh2const

## neutrino mass stuff.
dm21sqr = 7.62e-5
dm31sqr = 0.5*(2.55+2.43)*10**-3
Smnumin = np.zeros(2,dtype='float')
Smnumin[0] = dm31sqr**0.5 + dm21sqr**0.5
Smnumin[1] = dm31sqr**0.5 + np.sqrt(dm21sqr + dm31sqr)
print 'min Sum mnu', Smnumin

## fundamental constant.
DH = 2997.92458
kBoltzmann = 8.6173324e-5
myc1 = 120./7./(np.pi)**4.

## stuff needed to compute critical density of universe.
Mpc = 3.08567758e22
G = 6.67384e-11   #  m**3 kg** -1 s **-2
Msun = 1.9891e30  # kg
## critical density = 3.0*1.e10/G*Mpc/Msun/8./np.pi
rhocrit = 3.0*1.e10/G*Mpc/Msun/8./np.pi

def printsetup():
  global setup
  print setup

#############################
## begin class Pk        ####
#############################

## copying from cosmoxi2d/linearpk.c/getpkspline.

def Wtilde(x):
  """
  Fourier Transform of top-hat window fxn.
  This is written as a scalar function, could be vectorized.
  """
  xsqr = x*x
  if(x < 0.1):
    return (1.-xsqr/10.+xsqr*xsqr/280.-xsqr*xsqr*xsqr/15120.)
  else:
    return(3.*(np.sin(x)-x*np.cos(x))/(x*xsqr))


class linPk:
  def __init__(self,pkfname=None,ki=None,pki=None,z=0.,pkamprescale=1.):
    """
    Either input a camb file to read in, or a list of k and pk.  We set up pk interpolator.
    """
    if pkfname is not None:
      self.k, self.pk = np.loadtxt(pkfname,unpack=True,usecols=[0,1])
      self.pkfname = pkfname
    else:
      self.k = ki
      self.pk = pki
      self.pkfname = None
    if self.k is None or self.pk is None:
      self = None
      return None
    if len(self.k) != len(self.pk):
      self = None
      return None

    if np.fabs(pkamprescale - 1.) > 2.0e-3:
      self.pk = self.pk*pkamprescale

    self.npk = len(self.k)
    self.lnk = np.log(self.k)
    self.lnpk = np.log(self.pk)
    #self.pkinterp = interp.interp1d(self.lnk,self.lnpk,kind='cubic')
    ## equivalent to this!! but now we can have derivatives.  it's cubic by default.
    self.pkinterp = interp.InterpolatedUnivariateSpline(self.lnk,self.lnpk)
    self.nhigh = self.pkinterp.derivatives(self.lnk[-1])[1]
    self.nlow = self.pkinterp.derivatives(self.lnk[0])[1]

  def eval(self,kvec):
    lnkvec = np.log(np.fabs(kvec))
    try:
      len(lnkvec)
    except:
      lnkvec = np.array([lnkvec])

    assert len(lnkvec) > 0
    pkvec = np.zeros(len(lnkvec))
    xlow = np.where(lnkvec <= self.lnk[0])[0]
    xhigh = np.where(lnkvec >= self.lnk[-1])[0]
    x = np.where((lnkvec > self.lnk[0]) & (lnkvec < self.lnk[-1]))
    pkvec[x] = np.exp(self.pkinterp(lnkvec[x]))
    pkvec[xlow] = np.exp(self.lnpk[0] + self.nlow * (lnkvec[xlow] - self.lnk[0]))
    pkvec[xhigh] = np.exp(self.lnpk[-1] + self.nhigh * (lnkvec[xhigh] - self.lnk[-1]))
    return pkvec

  def sigmaR(self,R,verbose=False):
    ## want this to work for a np.array or for a scalar value.
    try:
      ll = len(R)
      Rtmp = R
      scalarflag = 0
    except:
      Rtmp = np.array([R])
      scalarflag = 1
    result = np.zeros(len(Rtmp))
    for i in range(len(Rtmp)):
      (ii, iierr) = quad(lambda k: k**2*self.eval(k)*(Wtilde(k*Rtmp[i]))**2,1.0e-8,10.)
      (ii2, iierr2) = quad(lambda k: k**2*self.eval(k)*(Wtilde(k*Rtmp[i]))**2,10.,500.)
      if verbose == True:
        print ii, iierr, ii2, iierr2
      result[i] = np.sqrt((ii+ii2)/2./np.pi**2)
    if scalarflag == 0:
      return result
    else:
      return result[0]


#############################
## end class Pk          ####
#############################


#############################
## begin exact neutrinos ####
#############################

def fexact(y):
  xmax = 100. ## tested limits, I think this is fine.
  ytmp = y
  scalarflag = 0
  try:
    ll = len(y)
  except:
    ytmp = np.array([y])
    scalarflag = 1
  result = ytmp.copy()
  for i in range(len(ytmp)):
    (ii, iierr) = quad(lambda x: x**2*np.sqrt(x**2 + (ytmp[i])**2)/(np.exp(x)+1.),0,xmax)
    result[i] = ii*myc1
  if(scalarflag == 1):
    return result[0]
  else:
    return result

def fKomatsu(y):
  A = 0.3173
  p = 1.83
  return ((1.+(A*y)**p)**(1./p))  

def makeftable():
  print 'making ftable..'
  global Ibigy
  global yi, fyi
  ymin = 0.
  ymax = 100.
  dy = 0.01
  yi=np.arange(ymin,ymax,dy)
  fyi = fexact(yi)
  ofp = open('ftable.dat','w')
  for i in range(len(yi)):
    ofp.write('%e %e\n' % (yi[i], fyi[i]))
  ofp.close() 
  Ibigy, tmperr = quad(lambda x: x**2/(np.exp(x)+1.),0,100)
  Ibigy = Ibigy*myc1 ## at large y, integral is Ibigy*y

def readftable():
  global Ibigy
  global yi, fyi
  yi, fyi = np.loadtxt('ftable.dat',unpack=True,usecols=[0,1])
  Ibigy, tmperr = quad(lambda x: x**2/(np.exp(x)+1.),0,100)
  Ibigy = Ibigy*myc1 ## at large y, integral is Ibigy*y
  return yi, fyi

## function defined by Komatsu that specifies energy density of neutrinos.
def f(y):
  global Ibigy
  global finterp
  try:  ## vector result.
    ll = len(y)
  except:  ## scalar result.
    if(y < yi[-1]):
      return finterp(y)
    else:
      return y*Ibigy

  result = y.copy()
  result[:] = Ibigy*y ## large y limit.
  ## for smaller y, use interpolation.
  xx = np.where(y < yi[-1])
  result[xx] = finterp(y[xx])
  return result

## copying from http://arxiv.org/pdf/1003.5918.pdf
def getm2m3(mlight,NorI):
  global dm21sqr, dm31sqr
  if(NorI == 0):
    m2 = np.sqrt(dm21sqr + mlight**2)
    m3 = np.sqrt(dm31sqr + mlight**2)
  else:
    m2 = np.sqrt(dm31sqr + mlight**2)
    m3 = np.sqrt(dm21sqr + m2**2)

  return m2, m3

def solveSmnu(mlight,NorI):
  m2, m3 = getm2m3(mlight,NorI)
  return mlight + m2 + m3


def setnumasses(Smnu,NorI):
  """
  Input sum of neutrino masses and hierarchy.  Output individual masses.
  """
  global dm21sqr, dm31sqr, Smnumin
  ## planck paper cites this paper:
  ## http://arxiv.org/pdf/1205.4018v4.pdf
  assert Smnu >= Smnumin[NorI]
  ## findroot for nu masses.
#  print solveSmnu(0,NorI), solveSmnu(Smnu/3.,NorI), Smnu
  mlight = scipy.optimize.brentq(lambda x: solveSmnu(x,NorI) - Smnu, 0, Smnu/3.)
  m2, m3 = getm2m3(mlight,NorI)
  #print 'hi beth',NorI,mlight,m2,m3,Smnu,(mlight+m2+m3)/Smnu
  return mlight, m2, m3

def onuh2tom1single(onuh2val,ogh2=2.469e-5,Tcmb=2.7255):
  """
  All the mass goes into a single species
  """
  Tnu_eV = Tcmb*((4./11.)**(1./3.))*kBoltzmann
  return (onuh2val/(onuh2const*ogh2*3.046/3.) - 2.*f(0.))/(f(100.)/100.)*Tnu_eV

def onuh2tomdegenerate(onuh2val,ogh2=2.469e-5,Tcmb=2.7255):
    """
    Mass split equally between 3 species
    """
    Tnu_eV = Tcmb*((4./11.)**(1./3.))*kBoltzmann
    return (onuh2val/(onuh2const*ogh2*3.046))/(f(100.)/100.)*Tnu_eV



#############################
##  end exact neutrinos  ####
#############################

def cosmosetup():
  global setup
  global yi, fyi
  global finterp
  if(setup == 1):
    return 0
  try:
    readftable()
  except:
    makeftable()
  finterp = interp.interp1d(yi,fyi,kind='linear')
  setup = 1
  return 0
  
## Planck XVI: Fixsen 2009 measured Tcmb = 2.7255 \pm 0.0006 K.

class cosmo:
  def __init__(self,och2=0.11219,obh2=0.02207,ogh2=2.469e-5,Tcmb=2.7255, \
                    h=0.7,w0=-1.,wa=0.,omegaDE=0.726,forceflat=1,\
                    Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
                    pkfname=None,zpk=None,pkamprescale=1.,MFopt=0,mffname=None):
    """
    This will store cosmological parameters and compute cosmological quantities of interest.
    Assume 3 std mass eigenstates (mnu1, mnu2, mnu3) and a single mass eigenstate for relativistic
    energy density for the (Neff-3.046) remaining neutrino-like particles.
    If you want to specify the neutrino mass sum and hierarchy and solve or the other 2 masses, 
    use NorI = 0 or 1 and specify Smnu > Smnumin.
    Or specify onuh2 and we solve for the mass assuming 
    either a single massive neutrino (SorD = 0) or 3 degenerate neutrinos (SorD = 1).
    Added functionality for Pk and stuff you might compute from that.
    For now you must input a camb file, in the future, we should make this thing call CAMB.
    In the future we should code up the growth equation or allow multiple P(k) files at different z.
    MFopt = 0: do nothing.  Leave MF blank.
    MFopt = 1: compute MF.  If mffname is not None, write it to that file. 
    MFopt = 2: read in MF from mffname.
    """

    cosmosetup()
    ## Temp and energy densities.
    self.Tcmb = Tcmb
    self.Tnu  = Tcmb*((4./11.)**(1./3.))
    self.Tnu_eV = self.Tnu*kBoltzmann

    self.och2 = och2
    self.obh2 = obh2
    self.ogh2 = ogh2
    self.h = h

    ## neutrino stuff.
    self.Neff = Neff #3.046/3. for each of mnu states 1-3, rest in mnu4
    self.mnu1 = mnu1
    self.mnu2 = mnu2
    self.mnu3 = mnu3
    self.mnu4 = mnu4
    assert(self.mnu4 < 1.0e-5) ## my tests on base_nnu_meffsterile_planck_lowl_lowLike_highL_post_lensing failed, I must have implemented something wrong for the sterile neutrino.
    ## optionally set from a sum and a hierarchy.
    if(NorI == 0 or NorI == 1):
      if(Smnu > Smnumin[NorI]):
        self.mnu1, self.mnu2, self.mnu3 = setnumasses(Smnu,NorI)
    else:
      if onuh2val > 0.:  ## convert to one mass eigenstate.
        if(SorD == 0):
          self.mnu1 = onuh2tom1single(onuh2val,ogh2=ogh2,Tcmb=Tcmb)
          assert(self.mnu2 < 1.0e-5)
          assert(self.mnu3 < 1.0e-5)
        if(SorD == 1):
          mval = onuh2tomdegenerate(onuh2val,ogh2=ogh2,Tcmb=Tcmb)
          self.mnu1 = mval
          self.mnu2 = mval
          self.mnu3 = mval
        assert np.fabs(self.onuh2(1.)/onuh2val - 1.) < 1.e-3 ## make sure only small numerical errors in the conversion.

    self.w0 = w0
    self.wa = wa
    self.oDE = omegaDE
    self.oDEh2 = self.oDE*h**2
    self.okh2 = h**2 - (self.och2 + self.obh2 + self.ogh2 + self.onuh2(1.) + self.oDEh2)
    if(forceflat==1):
      ## reset oDEh2 so that okh2 is 0.
      self.oDEh2 = self.oDEh2 + self.okh2
      self.oDE = self.oDEh2/h**2
      self.okh2 = 0.
      assert(1. - (self.och2 + self.obh2 + self.ogh2 + self.onuh2(1.) + self.oDEh2)/h**2 < 2.0e-6)
    self.ok = self.okh2/h**2
    if(np.fabs(self.ok) < 2.0e-6):
      self.flat = 1
    else:
      self.flat = 0

    ##derived parameters.
    ## non-neutrino matter density.
    self.ocb = (self.och2 + self.obh2)/(self.h)**2
    self.oc = (self.och2)/(self.h)**2
    self.ob = (self.obh2)/(self.h)**2
    self.og = (self.ogh2)/(self.h)**2
    #print 'deez',self.och2, self.obh2, self.ogh2, onuh2(1.,self)
    #print 'ok:',self.okh2,self.okh2/h**2

    if pkfname is not None and zpk is not None:
      self.linPk = linPk(pkfname=pkfname,pkamprescale=pkamprescale)
      self.zpk = zpk

    #MFopt = 1: compute MF.  If mffname is not None, write it to that file. 
    if MFopt == 1:
      if(mffname is not None):
        self.fillMF(mffname=mffname,writeopt = True)
      else:
        self.fillMF() # use defaults.

    #MFopt = 2: read in MF from mffname.
    if MFopt == 2:
      self.fillMF(mffname=mffname)


  def __str__(self):
    mystr = ''
    if(self.flat == 1):
      mystr += 'geometry: flat\n'
    else:
      mystr += 'geometry: omegak = %.3e\n' % self.ok
    mystr += 'h: %.4f\n' % self.h
    mystr += 'w0: %.4f\n' % self.w0
    mystr += 'wa: %.4f\n' % self.wa
    mystr += '\n'
    mystr += 'physical densities\n'
    mystr += 'och2: %.4f\n' % self.och2
    mystr += 'obh2: %.4f\n' % self.obh2
    mystr += 'ogh2: %.4e\n' % self.ogh2
    mystr += 'onuh2: %.4e\n' % self.onuh2(1.)
    mystr += 'oDEh2: %.4f\n' % self.oDEh2
    mystr += 'okh2: %.3e\n' % self.okh2
    mystr += '\n'
    mystr += 'neutrinos\n' 
    mystr += 'Neff: %.3e\n' % self.Neff
    mystr += 'masses: %.4f %.4f %.4f %.4f\n' % (self.mnu1, self.mnu2, self.mnu3, self.mnu4)
    return mystr

  ## this is meant to scale with a**-4 in H(a) = ...
  def onuh2(self,a):
    ## what to put here?
    ## all three masses.
    xf1 = self.mnu1*a/self.Tnu_eV
    xf2 = self.mnu2*a/self.Tnu_eV
    xf3 = self.mnu3*a/self.Tnu_eV
    xf4 = self.mnu4*a/self.Tnu_eV
    return onuh2const*self.ogh2*(3.046/3.*(f(xf1) + f(xf2) + f(xf3)) + \
                              (self.Neff-3.046)*f(xf4))
  def onu(self,a):
    ## what to put here?
    ## all three masses.
    xf1 = self.mnu1*a/self.Tnu_eV
    xf2 = self.mnu2*a/self.Tnu_eV
    xf3 = self.mnu3*a/self.Tnu_eV
    xf4 = self.mnu4*a/self.Tnu_eV
    return onuh2const*self.og*(3.046/3.*(f(xf1) + f(xf2) + f(xf3)) + \
                              (self.Neff-3.046)*f(xf4))

### need to generalize to non-flat, w0-wa etcetera models!!!  This is just to get going.
  def normH(self,a):
    """
    H(z)/(100 h km/s/Mpc)
    """
    ## copy from cosmoMCdigest
    myH2 = (self.ocb)*a**-3. + \
          self.oDE*a**(-3.*(1.+self.w0+self.wa))*np.exp(-3.*self.wa*(1.-a)) + \
          self.ok*a**-2. + \
          (self.og + self.onu(a))*a**-4.
    return np.sqrt(myH2)

  def normHinv(self,a):
    return(1./self.normH(a))

  def Hofz(self,a):
    return 100.*self.h*self.normH(a)

  def DofaLCDM(self,a):
    assert np.fabs(self.w0 + 1) < 2.0e-6
    assert np.fabs(self.wa) < 2.0e-6
    ## hmm, what about neutrinos?  make sure close to minimal at least.
    assert np.fabs(self.onuh2(1.) < 0.0007)
    i1, i1err = quad(lambda ap: 1./(ap*self.normH(ap))**3,1.0e-6,a)
    if np.fabs(i1err/i1) > 0.001:
      print 'growth convergence problem??'
      print i1, i1err
    ## not allowing neutrinos to contribute to omegam of growth.  correct?

    return (5./2.*self.ocb*self.normH(a)*i1)

  def dlnDdlnaLCDM(self,a):
    assert np.fabs(self.w0 + 1) < 2.0e-6
    assert np.fabs(self.wa) < 2.0e-6
    ## hmm, what about neutrinos?  make sure close to minimal at least.
    assert np.fabs(self.onuh2(1.) < 0.0007)
    
    return (1./(self.normH(a))**2*(-3./2.*self.ocb/a**3-self.ok/a**2 + 5./2.*a/self.DofaLCDM(a)*self.ocb/a**3))

  def dDdlnaLCDM(self,a):
    assert np.fabs(self.w0 + 1) < 2.0e-6
    assert np.fabs(self.wa) < 2.0e-6
    ## hmm, what about neutrinos?  make sure close to minimal at least.
    assert np.fabs(self.onuh2(1.) < 0.0007)

    myDofa = self.DofaLCDM(a)
    return myDofa/self.DofaLCDM(1.) * (1./(self.normH(a))**2*(-3./2.*self.ocb/a**3-self.ok/a**2 + 5./2.*a/myDofa*self.ocb/a**3))

  def DAz(self,a):
    """ this is in Mpc"""
    # converter to angular diameter distance.
    i1, i1err = quad(lambda ap: ap**-2*self.normHinv(ap),a,1)
    if(self.flat == 1):
      #print self.flat, self.ok, (i1*DH/self.h), (i1*DH/self.h) 
      return (i1*DH/self.h)
    else:  ## curvature
      if(self.ok > 0):
        #print self.flat, self.ok, (i1*DH/self.h), (DH/self.h/np.sqrt(self.ok)*np.sinh(np.sqrt(self.ok)*i1))
        return (DH/self.h/np.sqrt(self.ok)*np.sinh(np.sqrt(self.ok)*i1))
      else:
        #print self.flat, self.ok, (i1*DH/self.h), (DH/self.h/np.sqrt(-self.ok)*np.sin(np.sqrt(-self.ok)*i1)) 
        return (DH/self.h/np.sqrt(-self.ok)*np.sin(np.sqrt(-self.ok)*i1))

  def DVMpc(self,a):
    z = 1./a-1.
    return ((self.DAz(a))**2*self.normHinv(a)/self.h*DH*z)**(1./3.)

  def cambini(self,zoutlist,inifname,outbase):
    ofp = open(inifname,'w')
    ofp.write('DEFAULT(bethdefaults.ini)\n')
    ofp.write('output_root = %s\n' % (outbase))
    ofp.write('hubble = %.6f\n' % (self.h*100.))
    ofp.write('massive_neutrinos = 3.046\n')
    if(self.Neff > 3.046):
      ofp.write('massless_neutrinos=%.6f\n' % (self.Neff - 3.046))
    else:
      ofp.write('massless_neutrinos = 0.0\n')
    ofp.write('nu_mass_eigenstates = 3\n')
    ofp.write('nu_mass_degeneracies=1 1 1\n')
    mtot = (self.mnu1 + self.mnu2 + self.mnu3)
    if(mtot < 1.0e-5):
      mtot = 1.0e-5
      m1frac = 1./3.
      m2frac = 1./3.
      m3frac = 1./3.
    else:
      m1frac = self.mnu1/(self.mnu1 + self.mnu2 + self.mnu3)
      m2frac = self.mnu2/(self.mnu1 + self.mnu2 + self.mnu3)
      m3frac = self.mnu3/(self.mnu1 + self.mnu2 + self.mnu3)
    ofp.write('nu_mass_fractions = %.6e %.6e %.6e\n' % (m1frac, m2frac, m3frac))
    ofp.write('ombh2 = %.6f\n' % (self.obh2))
    ofp.write('omch2 = %.6f\n' % (self.och2))
    ofp.write('omnuh2 = %.6f\n' % (self.onuh2(1.)))
    ofp.write('transfer_num_redshifts  = %d\n' % (len(zoutlist)))
    for i in range(len(zoutlist)):
      ofp.write('transfer_redshift(%d) = %.6f\n' % (i+1,zoutlist[i]))
      ofp.write('transfer_filename(%d) = transfer_%d.out\n' % (i+1,i+1))
      ofp.write('transfer_matterpower(%d) = matterpower_%d.out\n' % (i+1,i+1))
    ofp.close()

##################################
## begin halo model utilities ####
## (mass fxn, bias fxn, ...)  ####
##################################

  def LagrangianMtoR(self,M):
    """
    Counts cdm and baryons towards mass budget, not massive neutrinos (!!??)
    Units are Msun/h.
    """
    denom = 4./3.*np.pi*rhocrit*self.ocb
    return (M/denom)**(1./3.)

  def LagrangianRtoM(self,R):
    """
    Comoving radius associated with 
    """
    return 4./3.*np.pi*rhocrit*self.ocb*R**3.

  def MtoRvir(self,MDelta,Delta=200.):
    return self.LagrangianMtoR(MDelta)*Delta**(-1./3.)

  def RtoMvir(self,RDelta,Delta=200.):
    return self.LagrangianRtoM(RDelta)*Delta

  def nuofM(self,M,deltac=1.686):
    return deltac/self.linPk.sigmaR(self.LagrangianMtoR(M))

  def nuofsigma(self,sigma,deltac=1.686):
    """
    Input sigma instead of M if you already computed sigma(M).
    """
    return deltac/sigma

  def bJTnu(self,nu,Delta=200.):
    y=np.log10(Delta)
    deltac = 1.686
    A = 1.0+0.24*y*np.exp(-(4./y)**4.)
    a = 0.44*y-0.88
    B = 0.183
    b = 1.5
    C = 0.019 + 0.107*y + 0.19*np.exp(-(4./y)**4.)
    c = 2.4
    return (1. - A*nu**a/(nu**a + deltac**a) + B*nu**b + C*nu**c)

  def bofM(self,M,Delta=200.):
    return self.bJTnu(nu=self.nuofM(M,Delta=Delta),Delta=Delta)

  def dndlnM(self,M,z,Delta=200.,wantsigma=False):
    """
    mass fxn at redshift of power spectrum you input.
    """
    if np.fabs(Delta - 200.) > 0.01:
      print 'I have not coded the interpolation for nofM; only currently works for Delta = 200.'
      return 0.
    if np.fabs(z - self.zpk) > 0.001:
      print 'code not generalized yet, mass fxn z must be hte same as power spectrum z.'
      return 0.
    ## parameters from http://adsabs.harvard.edu/abs/2008ApJ...688..709T for Delta = 200 (Table 2)
    A = 0.186*(1.+z)**-0.14
    a = 1.47*(1.+z)**-0.06
    alpha = 10**(-(0.75/np.log10(Delta/75.))**1.2)
    b = 2.57*(1.+z)**-alpha
    c = 1.19
    sig = self.linPk.sigmaR(self.LagrangianMtoR(M))
    ## make an interpolator object to get the log derivative.
    try:
      ll = len(M)
      ## ok, M is an array so proceed to derive derivative straight from there.
      lnsiginterp = interp.InterpolatedUnivariateSpline(np.log(M),np.log(sig))
      sigfac = np.array([-lnsiginterp.derivatives(np.log(mi))[1] for mi in M])
    except:
      ## M is a scalar, so evaluate sig at neighboring mass points and derive derivative.
      tmpfac = 0.9 ## factor (smaller than 1) to scale mass on either side to compute sigma derivative.
      mtmp = np.array([tmpfac*M,M,M/tmpfac])
      slist = np.array([self.linPk.sigmaR(self.LagrangianMtoR(tmpfac*M)), sig, self.linPk.sigmaR(self.LagrangianMtoR(M/tmpfac))])
      lnsiginterp = interp.InterpolatedUnivariateSpline(np.log(mtmp),np.log(slist))
      sigfac = -lnsiginterp.derivatives(np.log(M))[1]

    fsigma = A*((sig/b)**-a + 1.)*np.exp(-c/sig**2)
    if(wantsigma == False):
      return fsigma*rhocrit*self.ocb/M*sigfac
    else: ## return sigma vector corresponding to M as well.
      return sig, fsigma*rhocrit*self.ocb/M*sigfac

  def printMF(self,outfname):
    """
    Dump mass function theory to a file.
    """
    try:
      ll = len(self.mlist)
    except:
      return 0.

    ofp = open(outfname,'w')
    for i in range(len(self.mlist)):
      ofp.write('%e %e %e %e\n' % (self.mlist[i], self.siglist[i], self.mf[i], self.blist[i]))
    ofp.close()

  def fillMF(self,mmin=2.e11,mmax=1.e16,dlg10m=0.01,mffname=None,writeopt=False):
    """
    Evaluate the theory prediction for the mass function
    and compute/store sigma(M) along the way.
    Need to track Delta more carefully if it becomes a real variable.
    """
    if mffname is None or writeopt == True:
      ## log spacing.
      self.mlist = 10**(np.arange(np.log10(mmin),np.log10(mmax),dlg10m))
      self.siglist, self.mf = self.dndlnM(self.mlist,self.zpk, wantsigma=True)
      self.blist = self.bJTnu(nu=self.nuofsigma(self.siglist))  ## implicitly assumes Delta = 200. !!  carefull!!
    if writeopt == False and mffname is not None:
      self.mlist, self.siglist, self.mf, self.blist = np.loadtxt(mffname,unpack=True,usecols=[0,1,2,3])
    if writeopt == True and mffname is not None:
      self.printMF(mffname)

    self.mfinterp = interp.InterpolatedUnivariateSpline(np.log(self.mlist),np.log(self.mf))
    self.binterp = interp.InterpolatedUnivariateSpline(np.log(self.mlist),self.blist)

  def dndlnMinterp(self,M):
    return np.exp(self.mfinterp(np.log(M)))

  def biasinterp(self,M):
    return self.binterp(np.log(M))

##################################
## end halo model utilities   ####
## (mass fxn, bias fxn, ...)  ####
##################################



def Planck2cosmo(pelt,pnames):
  ## optional params
  ## fill with defaults.
  ## default for Planck is one mass eigenstate with mass 0.06: num_massive_neutrinos = 1
  ## however, if sum mnu is varying, then assume degenerate.
  popt = {'w':-1., 'wa':0.,'mnu':0.0588,'nnu':3.046,'meffsterile':0.,'omegak':0}
  for k in popt:
    if k in pnames:
      popt[k] = pelt[k]

  if 'omegak' in pnames:
    forceflat = 0
  else:
    forceflat = 1

  h = pelt['H0']*0.01
  #print 'this is h',h
  ocb = (pelt['omegach2'] + pelt['omegabh2'])/h**2
  onu = pelt['omeganuh2']/h**2
  #print pelt['omegamh2'], (pelt['omegach2'] + pelt['omegabh2']), pelt['omegamh2'] - (pelt['omegach2'] + pelt['omegabh2'] + pelt['omeganuh2'])
  ## what's in omegamh2??
  assert np.fabs(pelt['omegamh2'] - (pelt['omegach2'] + pelt['omegabh2'] + pelt['omeganuh2'])) < 2.0e-6
  #print ocb, pelt['omegal'],popt['omegak'],onu,ocb + pelt['omegal'] + popt['omegak'] + onu - 1.
  assert np.fabs(ocb + pelt['omegal'] + popt['omegak'] + onu - 1.) < 2.0e-6  ## this may fail, since radiation is 2e-5 ?

  ## to check: h, omegaDE = omegaL
  ## make sure for flat models I can get omega's to sum to 1
  ## makes sure I can derive onuh2 from mass inputs.

  ## two cases: minimal masses, all mass in single eigenstate.
  ## varying masses, assume degenerate.

  if('mnu' in pnames):
    cc = cosmo(och2=pelt['omegach2'],\
             obh2=pelt['omegabh2'],\
             h=pelt['H0']*0.01,\
             w0=popt['w'],\
             wa=popt['wa'],\
             omegaDE=pelt['omegal'],\
             forceflat = forceflat,\
             Neff=popt['nnu'],\
             mnu1=popt['mnu']/3.,\
             mnu2=popt['mnu']/3.,\
             mnu3=popt['mnu']/3.,\
             mnu4=popt['meffsterile'])
  else:
    ## minimal case.
    cc = cosmo(och2=pelt['omegach2'],\
             obh2=pelt['omegabh2'],\
             h=pelt['H0']*0.01,\
             w0=popt['w'],\
             wa=popt['wa'],\
             omegaDE=pelt['omegal'],\
             forceflat = forceflat,\
             Neff=popt['nnu'],\
             mnu1=popt['mnu'],\
             mnu2=0.,\
             mnu3=0.,\
             mnu4=popt['meffsterile'])

  ## sanity check on onuh2 now!
  if np.fabs(cc.onuh2(1.)/pelt['omeganuh2'] - 1.) > 1.0e-2:
    print 'neutrino mismatch??'
    print pelt['omeganuh2'], cc.onuh2(1.), onuh2tom1single(pelt['omeganuh2'],cc.ogh2,cc.Tcmb), onuh2tomdegenerate(pelt['omeganuh2'],cc.ogh2,cc.Tcmb), cc.mnu1, cc.mnu2, cc.mnu3

  #assert np.fabs(cc.onuh2(1.)/pelt['omeganuh2'] - 1.) < 1.0e-2
  return cc

def readPlanckchain(fparams, fchain):
  ## read in the parameter names.
  names = ['weight','lnlike']
  ifpp = open(fparams,'r')
  for line in ifpp:
    nn = line.split('\t')[0].strip(' ').strip('*')
    names.append(nn)

  chain = np.genfromtxt(fchain,names=names)
  return names,chain

def DAHcheckPlanck(fparams,fchain,nsamples=-1):
  """
  This function's purpose is to check my calculations 
  of DA and H at z=0.57 and DA(zstar) against the 
  values in the Planck chain.
  """
  names,chain = readPlanckchain(fparams,fchain)
  aval = 1./1.57
  maxdiff = 0.
  avgdiff = 0.
  avgDA = 0.
  maxdiffH = 0.
  avgdiffH = 0.
  avgH = 0.
  maxdiffDAstar = 0.
  avgdiffDAstar = 0.
  avgDAstar = 0.

  if nsamples == -1:
    nsamples = len(chain['omegach2'][:])

  myfac = len(chain['omegach2'][:])/nsamples

  i = 0
  while i < len(chain['omegach2'][:]):
#  for i in range(len(chain['omegach2'][:])):
    #if(i>10): break
    astar = 1./(1.+chain['zstar'][i])
    rstar = chain['rstar'][i]
    thetastar = chain['thetastar'][i]/100.
    DAstar = rstar/thetastar

    pcc = Planck2cosmo(chain[i],names)
    ## comparison between cc4 and myDA in LCDM demonstrates that at low redshift neutrinos can be treated like matter in terms of the expansion rate.
    ##cc4 = cosmo(och2=chain['omegach2'][i] + chain['omeganuh2'][i],obh2=chain['omegabh2'][i],h=chain['H0'][i]*0.01,forceflat=1)
    myDA = pcc.DAz(aval)*aval ## physical angular diameter distance.
    myDAstar = pcc.DAz(astar)  ## comoving angular diameter distance.  ARGH!
    myHofz = pcc.Hofz(aval)
    ##print myDA,myDA*aval,myDA*aval/chain['DA057'][i]-1,cc4.DAz(aval)/myDA-1,cc4.DAz(aval)*aval/chain['DA057'][i]-1
    #print myDA*aval,chain['DA057'][i],myDA*aval/chain['DA057'][i]-1,myDAstar, DAstar, myDAstar/DAstar-1.
    diff = np.fabs(myDA - chain['DA057'][i])
    avgdiff += diff
    maxdiff = max(maxdiff,diff)
    avgDA += myDA
    diffH = np.fabs(myHofz/DH/100. - chain['H057'][i])
    avgdiffH += diffH
    maxdiffH = max(maxdiffH,diffH)
    avgH += myHofz/DH/100.

    diffDAstar = np.fabs(myDAstar - DAstar)
    avgdiffDAstar += diffDAstar
    maxdiffDAstar = max(maxdiffDAstar,diffDAstar)
    avgDAstar += myDAstar
    i += myfac

  ll = float(len(chain['omegach2'][:]))
  avgdiff = avgdiff/ll
  avgdiffH = avgdiffH/ll
  avgdiffDAstar = avgdiffDAstar/ll
  avgDA = avgDA/ll
  avgDAstar = avgDAstar/ll
  avgH = avgH/ll
  print fchain.split('/')[-1],nsamples,ll,maxdiff/avgDA,maxdiffH/avgH, maxdiffDAstar/avgDAstar, \
                                 avgdiff/avgDA, avgdiffH/avgH, avgdiffDAstar/avgDAstar


def mysolveDAflat(h,cozold, astar, \
                   Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
                   w0=-1.,wa=0.):
  coztmp = cosmo(och2=cozold.och2, obh2=cozold.obh2, ogh2=cozold.ogh2, Tcmb=cozold.Tcmb, \
                 h=h,
                 w0=w0,
                 wa=wa,
                 forceflat=1,
                 Neff=Neff,mnu1=mnu1,mnu2=mnu2,mnu3=mnu3,mnu4=mnu4,NorI=NorI,
                 Smnu=Smnu,onuh2val=onuh2val,SorD=SorD)
  return coztmp.DAz(astar)
  #mlight = scipy.optimize.brentq(lambda x: solveSmnu(x,NorI) - Smnu, 0, Smnu/3.)

## this can be replaced with the curved version now, I checked they give the same answer.  I should keep the function name available though for backwards compatibility.
def cosmofromDAflat(cozold, DAzstar, zstar, \
                   Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
                   w0=-1.,wa=0.):
  """
  Input zstar.  Not allowed to change parameters that determine zstar; only vary H0.
  Enforces flatness.
  The purpose is to look at degeneracies at fixed DAstar.  
  Allowed to vary w0/wa or neutrino masses here from the fiducial model, and 
  h will be adjusted to keep DAzstar fixed.
  """
  astar = 1./(1.+zstar)
  try:
    newhval = scipy.optimize.brentq(lambda x: mysolveDAflat(x,cozold,astar,Neff, mnu1, mnu2, mnu3, mnu4, NorI, Smnu, onuh2val,SorD,w0,wa) - DAzstar,cozold.h*0.5,cozold.h*2.)
  ## should test for failure here!
  except:
    return None

  ## success!
  coznew = cosmo(och2=cozold.och2, obh2=cozold.obh2, ogh2=cozold.ogh2, Tcmb=cozold.Tcmb, \
                 h=newhval,
                 w0=w0,
                 wa=wa,
                 forceflat=1,
                 Neff=Neff,mnu1=mnu1,mnu2=mnu2,mnu3=mnu3,mnu4=mnu4,NorI=NorI,
                 Smnu=Smnu,onuh2val=onuh2val,SorD=SorD)
  assert((coznew.DAz(astar) - DAzstar)/DAzstar < 1.0e-5)
  return coznew

def mysolveDAcurved(h,cozold, astar, \
                   Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
                   w0=-1.,wa=0.,omegak=0.):

  if np.fabs(omegak) < 2.0e-6:
    forceflat = 1
  else:
    forceflat = 0
  okh2 = omegak*h**2
  oDEh2 = h**2 - (cozold.och2 + cozold.obh2 + cozold.ogh2 + cozold.onuh2(1.) + okh2)
  oDE = oDEh2/h**2
  coztmp = cosmo(och2=cozold.och2, obh2=cozold.obh2, ogh2=cozold.ogh2, Tcmb=cozold.Tcmb, \
                 h=h,
                 w0=w0,
                 wa=wa,
                 omegaDE = oDE,
                 forceflat=forceflat,
                 Neff=Neff,mnu1=mnu1,mnu2=mnu2,mnu3=mnu3,mnu4=mnu4,NorI=NorI,
                 Smnu=Smnu,onuh2val=onuh2val,SorD=SorD)
  return coztmp.DAz(astar)
  #mlight = scipy.optimize.brentq(lambda x: solveSmnu(x,NorI) - Smnu, 0, Smnu/3.)


def cosmofromDAcurved(cozold, DAzstar, zstar, \
                   Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
                   w0=-1.,wa=0.,omegak=0.):
  """
  Input zstar.  Not allowed to change parameters that determine zstar; only vary H0.
  Enforces flatness.
  The purpose is to look at degeneracies at fixed DAstar.  
  Allowed to vary w0/wa or neutrino masses here from the fiducial model, and 
  h will be adjusted to keep DAzstar fixed.
  """
  astar = 1./(1.+zstar)
  try:
    newhval = scipy.optimize.brentq(lambda x: mysolveDAcurved(x,cozold,astar,Neff, mnu1, mnu2, mnu3, mnu4, NorI, Smnu, onuh2val,SorD,w0,wa,omegak) - DAzstar,cozold.h*0.5,cozold.h*2.)
  ## should test for failure here!
  except:
    return None

  ## success!
  if np.fabs(omegak) < 2.0e-6:
    forceflat = 1
  else:
    forceflat = 0
  okh2 = omegak*newhval**2
  oDEh2 = newhval**2 - (cozold.och2 + cozold.obh2 + cozold.ogh2 + cozold.onuh2(1.) + okh2)
  oDE = oDEh2/newhval**2


  coznew = cosmo(och2=cozold.och2, obh2=cozold.obh2, ogh2=cozold.ogh2, Tcmb=cozold.Tcmb, \
                 h=newhval,
                 w0=w0,
                 wa=wa,
                 omegaDE = oDE,
                 forceflat=forceflat,
                 Neff=Neff,mnu1=mnu1,mnu2=mnu2,mnu3=mnu3,mnu4=mnu4,NorI=NorI,
                 Smnu=Smnu,onuh2val=onuh2val,SorD=SorD)
  assert((coznew.DAz(astar) - DAzstar)/DAzstar < 1.0e-5)
  return coznew


## same sort of excercise, but for w0-wa, fixing DAzstar and DAzeff.
#def 

## ran this over a Planck chain in my camb directory:
#howdiedoo@montserrat:~/mysoftware/cambMar13/camb$ python runfs8test.py > OUTrunfs8testplancklcdm

## Planck predicts fs8(z=0.57) = 0.481 \pm 0.010

def Planck2fs8camb(pelt, pnames, zhigh = 0.58, zlow = 0.56):
  if 'mnu' in pnames:
    print 'assuming here minimal neutrino mass, one mass eigenstate.'
    sys.exit(1)
  if zlow > zhigh:
    return -1

  cfname = 'cambfs8tmp.ini'
  cc = Planck2cosmo(pelt,pnames)
  print cc.DofaLCDM(1./1.56), cc.DofaLCDM(1./1.58), cc.dDdlnaLCDM(1./1.57)*pelt['sigma8'], pelt['sigma8']
  cfp = open(cfname,'w')
  cfp.write('ombh2 = %e\n' % cc.obh2)
  cfp.write('omch2 = %e\n' % cc.och2)
  cfp.write('omnuh2 = %e\n' % cc.onuh2(1.))
  cfp.write('hubble = %e\n' % (cc.h * 100))
  cfp.write('scalar_amp(1) = %e\n' % (np.exp(pelt['logA'])*1.e-10))
  cfp.write('scalar_spectral_index(1) = %e\n' % (pelt['ns']))
  cfp.write('DEFAULT(bethdefaultfs8.ini)\n')
  cfp.write('transfer_high_precision = T\n')
  cfp.write('transfer_kmax           = 2\n')
  cfp.write('transfer_k_per_logint   = 0\n')
  cfp.write('transfer_num_redshifts  = 3\n')
  cfp.write('transfer_interp_matterpower = T\n')
  cfp.write('transfer_redshift(3)    = 0.0\n')
  cfp.write('transfer_filename(3)    = transfer_z0p0.dat\n')
  cfp.write('transfer_redshift(2)    = %e\n' % zlow)
  cfp.write('transfer_filename(2)    = transfer_z0p56.dat\n')
  cfp.write('transfer_redshift(1)    = %e\n' % zhigh)
  cfp.write('transfer_filename(1)    = transfer_z0p58.dat\n')
  cfp.close() 
  ## now scrape the output.
  mystr = './camb %s > tmpout' % (cfname)
  os.system(mystr)

  sig8list = []
  ifp = open('tmpout','r')
  for line in ifp:
    if re.search('^ at z =',line):
      sig8val = float(line.split('(all matter)=')[-1].strip('\n'))
      sig8list.append(sig8val)
  ifp.close()
  if not len(sig8list) == 3:
    print sig8list
  assert len(sig8list) == 3
  ## check z=0
  if not np.fabs(sig8list[-1] - pelt['sigma8']) < 0.002:
    return -1 ## some alignment error.
  dlna = np.log((1.+zhigh)/(1.+zlow))
  fs8val = (sig8list[1] - sig8list[0])/dlna
  
  return fs8val
  




if __name__ == '__main__':
  printsetup()
  cosmosetup()
  printsetup()
  ## default cosmology used to measure the BOSS correlation function.
  cozdefault = cosmo()
  ## this looks correct, ./digestchainsbugfix/digestchainshackplanck code had 
  ## hi beth, fiducial values: 2.134761e+03 9.355810e+01 2.026613e+03 1.530572e+02
  ## which is good agreement (difference level 5e-5, which is probably because of here including radiation contribution)
  print 'default DA0p57:',cozdefault.DAz(1./1.57), cozdefault.DAz(1./1.57)/2.134761e+03-1.
  cc = cosmo(mnu1=0.1, mnu2=0.1, mnu3=0.1,forceflat=0)
  cc2 = cosmo(mnu1=0.1, mnu2=0.1, mnu3=0.1,forceflat=1)
  cc3 = cosmo(mnu1=0.0, mnu2=0.0, mnu3=0.0)
  cc4 = cosmo(Smnu = Smnumin[0]+0.0001,NorI = 0)

  zoutlist = np.array([4,2.,0])
  cc4.cambini(zoutlist,'cc4test.ini','cc4test')
  

  if(0==1): ## check the LCDM chains.
    pdir = '/home/howdiedoo/Planck/base_planck_lowl_lowLike_highL_post_lensing/base/planck_lowl_lowLike_highL/'
    pbase = 'base_planck_lowl_lowLike_highL_post_lensing'
    chainnum = 1
    DAHcheckPlanck(pdir+pbase+'.paramnames',pdir+pbase+'_%d.txt' % (chainnum))
  
  #fix ok problem. works now!  ok still the worst comparison, but it's ok.
  if(0==1): ## check the LCDM chains.
    pdir = '/home/howdiedoo/Planck/PLA/base_omegak/planck_lowl_lowLike_highL/'
    pbase = 'base_omegak_planck_lowl_lowLike_highL'
    chainnum = 1
    DAHcheckPlanck(pdir+pbase+'.paramnames',pdir+pbase+'_%d.txt' % (chainnum))

  if(0==1):
    pdir = '/home/howdiedoo/Planck/PLA/'
    pdir = '/Users/bareid/work/montserratdata/Planck/PLA/'
    chainnumlist = [1] ## just test one chain, not all of them.
    for chainnum in chainnumlist:
      os.system('ls %sbase*/planck_lowl_lowLike_highL/*highL_%d.txt > tmpoooq' % (pdir,chainnum))
      ifp = open('tmpoooq','r')
      for line in ifp:
        chainfname = line.strip('\n')
        mysplit = '_%d.txt' % (chainnum)
        pfname = chainfname.split(mysplit)[0] + '.paramnames'
        if(re.search('sterile',chainfname)):
          print 'skipping',chainfname,'for now.'
          continue
        DAHcheckPlanck(pfname,chainfname)
      ifp.close()

  if(0==1):  ## this one disagrees!  need to run this test and solve it, but will put this off til later.
    pdir = '/Users/bareid/work/montserratdata/Planck/PLA/base_nnu_meffsterile/planck_lowl_lowLike_highL/'
    pbase = 'base_nnu_meffsterile_planck_lowl_lowLike_highL_post_lensing' 
    chainnum = 1
    DAHcheckPlanck(pdir+pbase+'.paramnames',pdir+pbase+'_%d.txt' % (chainnum))

  ## works!!
  Smnu = 0.06
  while(Smnu < 1.2):
    if(Smnu > Smnumin[0]):
      m1,m2,m3 = setnumasses(Smnu,0)
      zeq = np.array([m1,m2,m3])/cc.Tnu_eV - 1.
#      print 'hi beth',0,'%.4f %.4f %.4f %.2f %.2f %.2f %.2f' % (m1,m2,m3,Smnu, zeq[0], zeq[1], zeq[2])
    if(Smnu > Smnumin[1]):
      m1,m2,m3 = setnumasses(Smnu,1)
      zeq = np.array([m1,m2,m3])/cc.Tnu_eV - 1.
#      print 'hi beth',1,'%.4f %.4f %.4f %.2f, %.2f %.2f %.2f' % (m1,m2,m3,Smnu, zeq[0], zeq[1], zeq[2])
    Smnu += 0.01

## good, we got Eqn 28 in Komatsu wmap_7yr_cosmology correct.
#  print 1./cc.Tnu_eV, 187.*3/94/10**-3, 187.*3/94/10**-3*cc.Tnu_eV
#  print cc.Tnu  ##also agreed.


