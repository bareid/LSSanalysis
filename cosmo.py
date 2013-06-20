import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

from scipy.integrate import quad
#import scipy.interpolate as interp

## fundamental constant.
DH = 2997.92458

### need to generalize to non-flat, w0-wa etcetera models!!!  This is just to get going.
def Eofz(z,om):
  return(np.sqrt(om*(1.+z)**3 + (1.-om)))

### om=0.274 is what we assume for all the mock catalogs and the dr9/dr10 datasets.
def DAz(z,om=0.274):
  # converter to angular diameter distance.
  return(quad(lambda zp: 1./Eofz(zp,om),0,z)[0]*DH)

class cosmo:
  def __init__(self,omegam=0.274,h=0.7):
    """
    This will store cosmological parameters and compute cosmological quantities of interest.
    """
    cosmo.om = omegam
    cosmo.h = 0.7


if __name__ == '__main__':
  rperpcut = 62./(180./np.pi*3600.)*DAz(0.7)
  print 'cutting rperp at ',rperpcut
  tmp = 62./(180./np.pi*3600.)
  print 'aaaa',tmp*DAz(0.43),tmp*DAz(0.55),tmp*DAz(0.7)


