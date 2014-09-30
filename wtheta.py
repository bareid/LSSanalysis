import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ximisc

class wtheta:
  def __init__(self,wtfname=None,icovfname=None,twlist=[],tmincut=-1.,tmaxcut=1.e12,skiprows=0):
    if((wtfname is not None) and (twlist == [])):
      self.fname = [wtfname] 
      try:
        t, wt = np.loadtxt(wtfname,skiprows=skiprows)
      except:
        print 'could not read wtfname'
        self = None
    else:
      try:
        t = twlist[0]
        wt = twlist[1]
        assert len(t) == len(wt)
      except:
        print 'bad format twlist'
        self = None

    if self is not None:
      self.t = t
      self.wt = wt
      self.ndata = len(self.t)
 
  ## finished init of wtheta!


def wtfromDR(fbase,fend='.ang'):
  """
  Quick hack, fill in options later.
  """

  fDD = fbase+'.DRopt1'+fend
  fDR = fbase+'.DRopt2'+fend
  fRR = fbase+'.DRopt3'+fend

  ximisc.getDRfactors(fbase,fend)

  t, DDg = np.loadtxt(fDD,skiprows=2,unpack=True)
  tchk, DRg = np.loadtxt(fDR,skiprows=2,unpack=True)
  assert (t==tchk).all()
  tchk, RRg = np.loadtxt(fRR,skiprows=2,unpack=True)
  assert (t==tchk).all()

  DRfac, fixRR = ximisc.getDRfactors(fbase,fend)

  wt = (DDg-DRg*DRfac)/RRg/DRfac**2/fixRR**2 + 1.

  return wtheta(wtfname=fbase,twlist = [t,wt])

def wtfromDDoDR(fbase,fend='.ang'):
  """
  Quick hack, fill in options later.
  """

  fDD = fbase+'.DRopt1'+fend
  fDR = fbase+'.DRopt2'+fend
  fRR = fbase+'.DRopt3'+fend

  DRfac = ximisc.getDRfactoronly(fbase,fend)

  t, DDg = np.loadtxt(fDD,skiprows=2,unpack=True)
  tchk, DRg = np.loadtxt(fDR,skiprows=2,unpack=True)
  assert (t==tchk).all()

  wt = 2.*DDg/DRg/DRfac - 1.

  return wtheta(wtfname=fbase,twlist= [t, wt])


