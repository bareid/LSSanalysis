import xiell
import wp
import xiwp
import numpy as np

class dset:
  def __init__(self,fbase,wdir,nbar2d=None,nbar3d=None,tasklist = [0,1,3],icovlist=[None,None,None]):
    """
    tasklist matches whichtask convention in catalog.callxi
    for each task, you can specify an icov matrix.
    """
    self.fbase = fbase
    self.wdir = wdir
    if 0 in tasklist:
      cut=0.533655
      ## remove the hard-coding here soon!
      try:
        indx = np.where(np.array(tasklist) == 0)[0]
        icovfname = icovlist[indx]
      except:
        icovfname = None
      self.xiell = xiell.xiellfromDR(wdir+"outputdr12-xiell/"+fbase,nell=2,binfile=wdir+"xibinfiles/bin1fineMU.txt",rperpcut=cut,icovfname=icovfname)
    if 1 in tasklist:
      try:
        indx = np.where(np.array(tasklist) == 1)[0]
        icovfname = icovlist[indx]
      except:
        icovfname = None
      self.wp = wp.wpfromDR(fbase=wdir+"outputdr12-xigrid/"+fbase,rpimax=80.0)
    if 3 in  tasklist:
      try:
        indx = np.where(np.array(tasklist) == 3)[0]
        icovfname = icovlist[indx]
      except:
        icovfname = None
      if nbar2d is None or nbar3d is None:
        print 'cannot compute wpcross without number densities!'
      else:
        self.nbar2d = nbar2d
        self.nbar3d = nbar3d
        self.wpcross = wp.wpcrossHogg(wdir+"outputdr12-wpcross/"+fbase,nbar2d=nbar2d,nbar3d=nbar3d)
        
        
def dsetNS(dsetN, dsetS, tasklist = [0,1,3], icovlist=[None,None,None]):
  x = dset(fbase=None,wdir=None,tasklist = []) ## do nothing.
  x.fbase = dsetN.fbase
  x.wdir = dsetN.wdir
  if 0 in tasklist:
    cut=0.533655
    ## remove the hard-coding here soon!
    try:
      indx = np.where(np.array(tasklist) == 0)[0]
      icovfname = icovlist[indx]
    except:
      icovfname = None
    tmp = dsetN.xiell + dsetS.xiell
    x.xiell = xiell.xiell(icovfname=icovfname, sxilist = [tmp.svec, tmp.xi])     
  if 1 in tasklist:
    try:
      indx = np.where(np.array(tasklist) == 1)[0]
      icovfname = icovlist[indx]
    except:
      icovfname = None
    tmp = dsetN.wp + dsetS.wp
    x.wp = wp.wp(icovfname=icovfname,rpwplist=[tmp.rsig, tmp.wp])
  if 3 in  tasklist:
    try:
      indx = np.where(np.array(tasklist) == 3)[0]
      icovfname = icovlist[indx]
    except:
      icovfname = None
    x.nbar2d = dsetN.nbar2d
    x.nbar3d = dsetN.nbar3d
    tmp = dsetN.wpcross + dsetS.wpcross
    x.wpcross = wp.wp(icovfname=icovfname,rpwplist=[tmp.rsig, tmp.wp])  
  return x     
