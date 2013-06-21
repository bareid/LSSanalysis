## this file contains miscellaneous useful functions.

import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

def downsample1dsum(a,dfac):
  adown = []
  for i in range(len(a)/dfac):
    aval = (a[i*dfac:(i+1)*dfac]).sum()
    adown.append(aval)
  adown = np.array(adown)
  return adown


def downsample1d(a,dfac):
  adown = []
  for i in range(len(a)/dfac):
    aval = (a[i*dfac:(i+1)*dfac]).sum()
    adown.append(aval)
  adown = np.array(adown)
  adown = adown/float(dfac)
  return adown

def legendre(n,x):
  if(n==0):
    return 1
  if(n==1):
    return x
  if(n==2):
    return 0.5*(3.*x**2-1.)
  if(n==3):
    return 0.5*(5.*x**3 -3.*x)
  if(n==4):
    return 0.125*(35.*x**4 -30.*x**2 + 3.)
  if(n>4 or n<0):
    sys.exit(1)

def one2twod(a,b):
  a2d = np.zeros([len(a), len(b)])
  b2d = np.zeros([len(a), len(b)])
  for i in range(len(a)):
    a2d[i,:] = a[i]
    b2d[i,:] = b

  return a2d, b2d

## generate comparison to see if two results are close (enough).
#these can be any form of xi/wp/xiell, as long as they have the same.
def comparexi(xia, xib,maxfrac=1.e-4,maxdiff=1.e-4):
  """
  maxdiff is the maximum absolute difference allowed.
  maxfrac is maximum fractional deviation tolerated (unless diff less than maxdiff)
  Returns 0 if they're the same.
  Returns 1 if they disagree.
  Returns 2 if they have different shapes.
  """
  xi1 = xia.flatten()
  xi2 = xib.flatten()
  if len(xi1) != len(xi2):
    return 2
  xx = np.where(np.fabs((xi1-xi2)/xi1) > maxfrac,1,0)
  yy = np.where(np.fabs(xi1-xi2) > maxdiff,1,0)
  zz = (xx & yy)

  if (zz != 0).any():
    return 1
  else:
    return 0


def comparegeneric(a, b,maxfrac=1.e-4,maxdiff=1.e-4):
  """
  maxdiff is the maximum absolute difference allowed.
  maxfrac is maximum fractional deviation tolerated (unless diff less than maxdiff)
  Returns 0 if they're the same.
  Returns 1 if they disagree.
  Returns 2 if they have different shapes.
  This is identical to comparexi, should I delete that one?
  """
  xi1 = a.flatten()
  xi2 = b.flatten()
  if len(xi1) != len(xi2):
    return 2
  xx = np.where(np.fabs((xi1-xi2)/xi1) > maxfrac,1,0)
  yy = np.where(np.fabs(xi1-xi2) > maxdiff,1,0)
  zz = (xx & yy)

  if (zz != 0).any():
    return 1
  else:
    return 0

def getmatrixnorm(m):
  mnorm = copy.deepcopy(m)
  (nx, ny) = m.shape
  for i in range(nx):
    for j in range(ny):
      mnorm[i,j] = m[i,j]/(m[i,i]*m[j,j])**0.5 
  return mnorm

def getmatrixdiag(m):
  (nx, ny) = m.shape
  if nx != ny:
    return None
  diag = np.zeros(nx)
  for i in range(nx):
    diag[i] = m[i,i]
  return diag

##read out D/R factors from pair counts headers.
def getDRfactors(fbase):
  """
  Input DR and RR filepaths.  This function reads
  DRfac and fixRRdown from that path.
  """

  fDR = fbase+'.DRopt2'
  fRR = fbase+'.DRopt3'

  ifpRR = open(fRR,'r')
  line = ifpRR.readline()
  line = ifpRR.readline()
  DDwgt = float(line.split(':')[1].split(',')[0])
  RRwgt = float(line.split(':')[1].split(',')[1])
  RRfacdown = float(RRwgt)
  ifpRR.close()

  ifpDR = open(fDR,'r')
  line = ifpDR.readline()
  line = ifpDR.readline()
  ## rerun everything to get this to 12 digits, just in caes!
  DDwgt = float(line.split(':')[1].split(',')[0])
  RRwgt = float(line.split(':')[1].split(',')[1])
  DRfac = float(DDwgt)/float(RRwgt)
  ifpDR.close()
  fixRR = RRwgt/RRfacdown

  return DRfac, fixRR

def getDRfactoronly(fbase):
  """
  Input DR filepaths.  This function reads
  DRfac and fixRRdown from that path.
  """

  fDR = fbase+'.DRopt2'

  ifpDR = open(fDR,'r')
  line = ifpDR.readline()
  line = ifpDR.readline()
  ## rerun everything to get this to 12 digits, just in caes!
  DDwgt = float(line.split(':')[1].split(',')[0])
  RRwgt = float(line.split(':')[1].split(',')[1])
  DRfac = float(DDwgt)/float(RRwgt)
  ifpDR.close()

  return DRfac
