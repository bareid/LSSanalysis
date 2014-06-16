import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import cosmo

class halocat:
  def __init__(self,fname=None,mindx=None,mlist=[],logopt = 0, Lbox=2750., massfxnfname=None):
    """
    Read in a list of masses.  If logopt==1, then input masses are understood to be logarithmic.
    """
    self.Lbox = Lbox
    if mlist != []:
      if logopt == 0:
        self.m = np.array(mlist)
        self.lg10m = np.log10(self.m)
      else:
        self.lg10m = np.array(mlist)
        self.m = 10**(self.log10m)
      self.lg10mcen, self.Nofm = None, None ## set these later with massfxn.
    elif massfxnfname is not None:
      self.lg10mcen, self.Nofm = np.loadtxt(massfxnfname,unpack=True,usecols=[0,1])


    else:
      if 0==0:
#      try:
        if logopt == 0:
          self.m = np.loadtxt(fname,usecols=[mindx],unpack=True)
          self.lg10m = np.log10(self.m)
        else:
          self.lg10m = np.loadtxt(fname,usecols=[mindx],unpack=True)
          self.m = 10**(self.lg10m)
        self.lg10mcen, self.Nofm = None, None ## set these later with massfxn.
      else:
#      except:
        print 'file read did not work.'
        self.m = None
        self.lg10m = None

  def massfxn(self,dlg10m=0.01,lg10mmin=None,lg10mmax=None):
    """ Compute (log binned) mass fxn """
    if(lg10mmin is None or lg10mmax is None):
      h,x = np.histogram(self.lg10m,bins=np.arange(self.lg10m.min()-0.5*dlg10m, self.lg10m.max()+0.55*dlg10m,dlg10m))
    else:
      h,x = np.histogram(self.lg10m,bins=np.arange(lg10mmin, lg10mmax, dlg10m))

    hnorm = h/dlg10m/self.Lbox**3
    self.lg10mcen = 0.5*(x[1:]+x[:-1])
    self.Nofm = hnorm
#    print 'yo',self.lg10mcen, self.Nofm

  def addmassfxncurve(self,ax=None,dlg10m=0.01,lg10mmin=None,lg10mmax=None,color='k',lbl=None):
    """
    Add mass fxn to a curve.
    """
    self.massfxn(dlg10m,lg10mmin,lg10mmax)
    if ax is None:
      plt.plot(10**(self.lg10mcen),self.Nofm,color=color,label=lbl)
    else:
      ax.plot(10**(self.lg10mcen),self.Nofm,color=color,label=lbl)

  def massfxnplot(self,dlg10m=0.01,lg10mmin=None,lg10mmax=None,color='k',lbl=None):
    """ Plot a mass fxn"""
    f=plt.figure()
    ax = f.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    self.addmassfxncurve(ax,dlg10m,lg10mmin,lg10mmax,color,lbl)
    return f, ax

class halo2pt:
  def __init__(self,fname=None):
    """
    Read in 2pt xi and v stats.  Compare cats.
    """
    try:
      r, xi, vr, v2par, v2perp = np.loadtxt(fname,unpack=True)
      skiprows = max((np.where(np.isnan(vr))[0]).max(), (np.where(np.isnan(v2par))[0]).max(), (np.where(np.isnan(v2perp))[0]).max()) + 1
      self.r, self.xir, self.vr, self.v2par, self.v2perp = np.loadtxt(fname,unpack=True,skiprows=skiprows)
      self.nr = len(self.r)
      ## comput binning.
      self.logopt = -1 ## non-simple binning.
      dr = (self.r[1:]-self.r[:-1]).mean()
      if (np.fabs(self.r[1:]-self.r[:-1] - dr) < 0.0001*dr).all():
        self.logopt = 0
        self.dr = dr
      dlogr = (np.log(self.r[1:]/self.r[:-1])).mean()
      if (np.fabs(np.log(self.r[1:])-np.log(self.r[:-1]) - dlogr) < 0.0001*dlogr).all():
        self.logopt = 1
        self.dlogr = dlogr


    except:
      self = None

  ## copying from wpinterp
  def xirinterp(self,rarr):
    ## do interpolation in log space.
    if ((type(rarr) is type(0.0)) | (type(rarr) is type(0))):
      rarr = np.array([rarr])

    ## make sure it's type numpy
    rarr = np.array(rarr)

    if(self.logopt == -1):
## need to write this piece later.
      return 1.

    if(self.logopt == 0):
      ix = (rarr-self.r.min())/self.dr 
    if(self.logopt == 1):
      ix = (np.log(rarr/self.r.min()))/self.dlogr 
    tmp = np.where(ix < 0)[0]
    ix[tmp] = 0
    tmp = np.where(ix > (self.nr)-1)[0]
    ix[tmp] = self.nr-1.001
    iix = np.array([int(val) for val in ix])
    assert (iix >= 0).all()
    assert (iix <= self.nr-2).all()
    fx = ix -iix
    assert (fx >= 0.).all()
    assert (fx <= 1.).all()


    if(self.logopt == 0):
      if(len(iix) > 1):
        return self.xir[iix]*(1.-fx) + self.xir[iix+1]*fx
      else:
        return (self.xir[iix]*(1.-fx) + self.xir[iix+1]*fx)[0]
    if(self.logopt == 1):
      if(len(iix) > 1):
        return np.exp(np.log(self.xir[iix])*(1.-fx) + np.log(self.xir[iix+1])*fx)
      else:
        return (np.exp(np.log(self.xir[iix])*(1.-fx) + np.log(self.xir[iix+1])*fx))[0]

def plothalo2pt(hlist,clist,fside,ci,normr=10.):
  """
  Plot xi(r), v(r) and dispersions.
  ci = central index, the one you want to make ratios/differences relative to.
  """

  figsize=[3.*fside,2.*fside]
  f = plt.figure(figsize=figsize)
  ax1 = f.add_subplot(231)
  ax2 = f.add_subplot(232)
  ax3 = f.add_subplot(233)
  ax4 = f.add_subplot(234)
  ax5 = f.add_subplot(235)
  ax6 = f.add_subplot(236)
  for c, h in zip(clist,hlist):
    nstart = 0
    dstart = 0
    if(len(hlist[ci].xir) > len(h.xir)):
      dstart = len(hlist[ci].xir) - len(h.xir)
    if(len(hlist[ci].xir) < len(h.xir)):
      nstart = -len(hlist[ci].xir) + len(h.xir)

    xinormfac = (h.xirinterp(normr)/hlist[ci].xirinterp(normr))
    vnormfac = xinormfac**0.5
    print 'normfacs for r = ',normr,':',xinormfac, vnormfac

#    print 'yo',nstart, dstart, len(hlist[ci].xir), len(h.xir)
    ax1.plot(h.r, h.xir,color=c)
    ax4.plot(h.r[nstart:], h.xir[nstart:]/hlist[ci].xir[dstart:]/xinormfac,color=c,linestyle='--')
    ax4.plot(h.r[nstart:], h.xir[nstart:]/hlist[ci].xir[dstart:],color=c)
    ax2.plot(h.r, h.vr,color=c)
    ax5.plot(h.r[nstart:], h.vr[nstart:]/hlist[ci].vr[dstart:],color=c)
    ax5.plot(h.r[nstart:], h.vr[nstart:]/hlist[ci].vr[dstart:]/vnormfac,color=c,linestyle='--')
    ax3.plot(h.r, (h.v2perp),color=c,ls='-')
    ax3.plot(h.r, (h.v2par),color=c,ls='--')
    diff = h.v2perp[-1] - hlist[ci].v2perp[-1]
    ax6.plot(h.r, (h.v2perp)-diff,color=c,ls='-')
    diff = h.v2par[-1] - hlist[ci].v2par[-1]
    ax6.plot(h.r, (h.v2par)-diff,color=c,ls='--')

  for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_xscale('log')

  return f, [ax1, ax2, ax3, ax4, ax5, ax6]

if __name__ == '__main__':

  ### match mass fxn done in precompute.c
  dlg10m = 0.01
  hhfof1 = halocat(fname="/home/howdiedoo/SOforL0work/packSOforL0openmp/SOconvertboundary/FOFnomatchfile.dat.cut",mindx=7,logopt=1)
  lg10mmin = hhfof1.lg10m.min()-0.5*dlg10m

  ## tmp!
  #hhfof1=halocat(fname="/home/howdiedoo/SOforL0work/packSOforL0openmpnewvel/SOforL0.concat.cut",mindx=7,logopt=1)
  #lg10mmin = 11.
  ## end tmp!!
  print 'min',lg10mmin
  hhfofp=halocat(fname="/home/howdiedoo/SOforL0work/L00_0.6452.halos",mindx=8,logopt=1)
#  hhfofp=halocat(fname="/home/howdiedoo/SOforL0work/packSOforL0openmpnewvel/SOforL0.concat.cut",mindx=7,logopt=1)
#  hhSO=halocat(fname="/home/howdiedoo/SOforL0work/packSOforL0openmpnewvel/SOforL0.concat.cut",mindx=7,logopt=1)
  hhSO=halocat(fname="/home/howdiedoo/SOforL0work/packSOforL0openmpnewvel/SOforL0.concat",mindx=7,logopt=1)
  lg10mmax = max(hhfof1.lg10m.max(),hhfofp.lg10m.max(),hhSO.lg10m.max())+0.51*dlg10m
  print 'max',lg10mmax
  hhfof1.massfxn(dlg10m,lg10mmin,lg10mmax)
  hhSO.massfxn(dlg10m,lg10mmin,lg10mmax)
  hhfofp.massfxn(dlg10m,lg10mmin,lg10mmax)

  print 'ga',hhfof1.lg10mcen
  assert (hhfof1.lg10mcen == hhSO.lg10mcen).all()
  assert (hhfof1.lg10mcen == hhfofp.lg10mcen).all()

  ## print out hte mass fxns.
  ofp = open("L0massfxns.dat",'w')
  for i in range(len(hhSO.lg10mcen)):
    ofp.write('%f %e %e %e\n' % (hhSO.lg10mcen[i], hhfofp.Nofm[i], hhfof1.Nofm[i], hhSO.Nofm[i]))
    
  ofp.close()

  f, ax = hhfofp.massfxnplot(dlg10m,lg10mmin,lg10mmax,color='k',lbl='FOFp')
  hhfof1.addmassfxncurve(ax,dlg10m,lg10mmin,lg10mmax,color='b',lbl='FOF1')
  hhSO.addmassfxncurve(ax,dlg10m,lg10mmin,lg10mmax,color='g',lbl='SO')
  ax.plot(hhSO.lg10mcen, hhSO.Nofm+hhfof1.Nofm,color='r',lbl='sum')
  plt.legend(loc=3)
  f.savefig("massfxnall3.png")
  


