import numpy as np
import os
import re

def diffch(dir1,dir2,outfile=None):
  """ appends diff outputs to outfile"""
  for ff in sorted(os.listdir(dir1)):
    if re.search('.c$',ff) or re.search('.h$',ff):
      f1 = dir1 + ff
      f2 = dir2 + ff
      if outfile is None:
        print 'start diff ',f1,f2
        os.system('diff %s %s' % (f1,f2))
        print 'end   diff ',f1,f2
      else:
        ofp = open(outfile,'a')
        ofp.write('start diff %s %s\n' % (f1,f2))
        ofp.close()
        os.system('diff %s %s >> %s' % (f1,f2,outfile))
        ofp = open(outfile,'a')
        ofp.write('end   diff %s %s\n' % (f1,f2))
        ofp.close()




def rebin(v, wgt = None,istart=0,rebinsize=3):
  """
  reduce the fineness of the binning by factor rebinsize, starting at bin istart.
  optionally, use a weight for each bin (rather than assuming equal weight in the averaging).
  """
  if wgt is None:
    wgt = np.zeros(len(v),dtype='float') + 1.

  vdown = np.zeros((len(v)-istart)/rebinsize)

  ii = istart
  iidown = 0
  while ii <= len(v) - rebinsize:
    vdown[iidown] = (wgt[ii:ii+rebinsize]*v[ii:ii+rebinsize]).sum()/(wgt[ii:ii+rebinsize]).sum()
    iidown += 1
    ii += rebinsize

  return vdown

def rebinedge(vx, vy, vxedge, wgt = None, nonorm=False):
  if wgt is None:
    wgt = np.zeros(len(vx),dtype='float') + 1.

  vxdown = np.zeros(len(vxedge)-1)
  vydown = np.zeros(len(vxedge)-1)
  wgtdown = np.zeros(len(vxedge)-1)

  for ii in range(len(vxdown)):
    xx = np.where((vx >= vxedge[ii]) & (vx < vxedge[ii+1]))[0]
    if nonorm == True:
      vxdown[ii] = (wgt[xx]*vx[xx]).sum()
      vydown[ii] = (wgt[xx]*vy[xx]).sum()
      wgtdown[ii] = (wgt[xx]).sum()

    else:
      vxdown[ii] = (wgt[xx]*vx[xx]).sum()/(wgt[xx]).sum()
      vydown[ii] = (wgt[xx]*vy[xx]).sum()/(wgt[xx]).sum()

  if nonorm == True:
    return vxdown, vydown, wgtdown
  else:
    return vxdown, vydown  ## for compatibility with old version.
