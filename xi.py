import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import re
import xi2d
import wp
import xiell


class xi:
  def __init__(self,xi2dfname=None,wpfname=None,wpicovfname=None,xiellfname=None,xiellicovfname=None,xidatfname=None):
    """
    Uses file inputs to read in various correlation function statistics (xi2d, xiell, wp) into a xi object.
    The main purpose of this class is to package a single model/measurement together and make a fancy plot.
    """

    if xidatfname is not None:
      ifp = open(xidatfname,'r')
      print 'reading info from',xidatfname
      for line in ifp:
        if xi2dfname is None:
          if(re.match('xi2d:',line)):
            xi2dfname = line.split(':')[1].strip(' \n')
        if wpfname is None:
          if(re.match('wp:',line)):
            x = line.split(':')[1].split(',')
            assert len(x) == 1 or len(x) == 2
            wpfname = x[0].strip(' \n')
            if(len(x) == 2 and wpicovfname is None):
              wpicovfname = x[1].strip(' \n')
        if xiellfname is None:
          if(re.match('xiell:',line)):
            x = line.split(':')[1].split(',')
            assert len(x) == 1 or len(x) == 2
            xiellfname = x[0].strip(' \n')
            if(len(x) == 2 and xiellicovfname is None):
              xiellicovfname = x[1].strip(' \n')

    self.xi2dfname = xi2dfname
    self.wpfname = wpfname
    self.xiellfname = xiellfname

    if(xi2dfname is not None):
      self.xi2d = xi2d.xi2d(xi2dfname)
    else: self.xi2d = None
    if(wpfname is not None):
      self.wp = wp.wp(wpfname,icovfname=wpicovfname)
    else: self.wp = None
    if(xiellfname is not None):
      self.xiell = xiell.xiell(xiellfname,icovfname=xiellicovfname)
    else: self.xiell = None

  def addfancyplot(self,axxi2d,axwp,axxiell,settings):
#color='k',fmt=None,lbl=None,\
#    clevlist=[],symmetrizeopt=1,\
#    wppow=0,xiellpow=1,elllist=[0,2],color2=None):
    """
    add xi2d, xiell, wp curves to a fancy plot.
    """
    if settings['color2'] is None:
      color2 = settings['color']
    else:
      color2 = settings['color2']

    if(self.xi2d is not None):
      self.xi2d.addcontour(ax=axxi2d,symmetrizeopt=settings['symmetrizeopt'],\
        clevlist=settings['clevlist'],color=settings['color'])
    if(self.wp is not None):
      self.wp.addcurve(ax=axwp,color=settings['color'],rppow=settings['wppow'],\
        fmt=settings['fmt'],lbl=settings['lbl'])
    if(self.xiell is not None):
      for ellval in settings['elllist']:
        if(ellval == 0):
          self.xiell.addcurve(ax=axxiell,ell=ellval,color=settings['color'],\
            spow=settings['xiellpow'],fmt=settings['fmt'],lbl=settings['lbl'])    
        else:
          self.xiell.addcurve(ax=axxiell,ell=ellval,color=color2,\
            spow=settings['xiellpow'],fmt=settings['fmt'],lbl=settings['lbl']) 

  def makefancyplot(self,sizescale=1.,sizeratio=2.,color='k',fmt=None,lbl=None,\
         clevlist=[],symmetrizeopt=1,spanxi2d=None,\
         wppow=0,spanwp=None,logxoptwp=1,logyoptwp=1,\
         xiellpow=1,spanxiell=None,logxoptxiell=1,logyoptxiell=0,elllist=[0,2],color2=None,customax=True):
        
    """
    make fancy plot.  aspect ratios are fixed, 
    but you can scale the total size of the resulting plot [sizescale]
    or the sizeratio between the big (xi2d) and small (wp, xiell) plots [sizeratio].
    color,color2[used for xi2], fmt,lbl are passed to all the plots.  
    the rest of the options are plot-specific and labelled that way.
    Returns the axes objects, fig object is an attribute of the class.
    """
    
    ## define some buffers on the left, right, top bottom, and in between.
    ## if you want aspect ratio = 1, these better sum up to the same values (xtot and ytot)!
    xbl = 0.1
    xbm = 0.1
    xbr = 0.03
    ## eqns below break if you don't set x and y borders equal
    ybt = xbr
    ybm = xbm
    ybb = xbl
    dy2d=0.5*ybm
    
    ## ratio of size of xi2d plot to wp/xiell plots.
    sizeratio = 2. 
    ## space remaining after whitespace.
    xtot = 1.-xbl-xbm-xbr
    ytot = 1.-ybt-ybm-ybb
    xsize2d = (xtot)*sizeratio/(sizeratio+1.)
    xsize1d = (xtot)*1./(sizeratio+1.)
    ysize2d = xsize2d
    ysize1d = xsize1d
    
    #all the y distances need to be scaled by 1/yfac
    yfac = (ybb + ybm + ybt + ysize2d)/(ybb + ybm + ybt + ysize2d+ysize1d)
    
    
#    print 'stuff',xsize2d, xsize1d, dy2d
#    print xbl+xsize2d, ybb+dy2d+ysize2d
#    print 'axxi2d: ',xbl,ybb+dy2d,xsize2d,ysize2d
    
    
    self.fancyfig = plt.figure(figsize=[sizescale*12,sizescale*12*yfac])
    axxi2d = self.fancyfig.add_axes([xbl,(ybb+dy2d)/yfac,xsize2d,ysize2d/yfac])
    
    x2 = xbl + xsize2d+xbm
    ywp = ybb + ysize1d+ybm
#    print 'wp: ',x2,ywp,xsize1d,ysize1d
    axwp = self.fancyfig.add_axes([x2,ywp/yfac,xsize1d,ysize1d/yfac])
#    print 'xiell: ',x2,ybb,xsize1d,ysize1d
#    print 'gaa',x2+xsize1d,ywp+ysize1d
    axxiell = self.fancyfig.add_axes([x2,ybb/yfac,xsize1d,ysize1d/yfac])

    if(logxoptwp==1): 
      axwp.set_xscale('log')
    # else: don't set to linear, just leave as it is.
    if(logyoptwp==1):
      axwp.set_yscale('log')
    if spanwp is None:
      spanwp = [self.wp.rsig.min()*0.9, self.wp.rsig.max()*1.1, (self.wp.wp*self.wp.rsig**wppow).min()*0.9, (self.wp.wp*self.wp.rsig**wppow).max()*1.1]
    axwp.axis(spanwp)
    
        
    if(logxoptxiell==1): 
      axxiell.set_xscale('log')
    # else: don't set to linear, just leave as it is.
    if(logyoptxiell==1):
      axxiell.set_yscale('log')    
    if spanxiell is None:
      ell=0
      spanxiell0 = [self.xiell.svec[ell/2,:].min()*0.9,self.xiell.svec[ell/2,:].max()*1.1,\
      (self.xiell.svec[ell/2,:]**xiellpow*self.xiell.xi[ell/2,:]).min()*0.9,\
      (self.xiell.svec[ell/2,:]**xiellpow*self.xiell.xi[ell/2,:]).max()*1.1]
      ell=2
      spanxiell2 = [self.xiell.svec[ell/2,:].min()*0.9,self.xiell.svec[ell/2,:].max()*1.1,\
      (self.xiell.svec[ell/2,:]**xiellpow*self.xiell.xi[ell/2,:]).min()*0.9,\
      (self.xiell.svec[ell/2,:]**xiellpow*self.xiell.xi[ell/2,:]).max()*1.1]
      spanxiell=copy.copy(spanxiell0)
      spanxiell[0] = min(spanxiell0[0], spanxiell2[0])
      spanxiell[2] = min(spanxiell0[2], spanxiell2[2])
      spanxiell[1] = max(spanxiell0[1], spanxiell2[1])
      spanxiell[3] = max(spanxiell0[3], spanxiell2[3])

    axxiell.axis(spanxiell)

    if spanxi2d is not None:
      axxi2d.axis(spanxi2d)
        
    fancydict= {'color':color,'fmt':fmt,'lbl':lbl,'clevlist':clevlist,\
     'symmetrizeopt':symmetrizeopt,'wppow':wppow,'xiellpow':xiellpow,\
     'elllist':elllist,'color2':color2}
    self.fancydict = copy.deepcopy(fancydict) ## remember the values!

    self.addfancyplot(axxi2d,axwp,axxiell,settings=fancydict)
#color,fmt,lbl,clevlist,symmetrizeopt,wppow,xiellpow,elllist,color2,fancyplotsettings)

    if(customax==True):
      customxticks=[0.5,1.0,5.0,10.,20.,30.]
      axwp.xaxis.set_major_locator(plt.FixedLocator(customxticks))
      axwp.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
      axwp.xaxis.set_ticks_position('bottom')
      axwp.tick_params(axis='x',reset=False,which='both',length=8,width=2)
      axxiell.xaxis.set_major_locator(plt.FixedLocator(customxticks))
      axxiell.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
      axxiell.xaxis.set_ticks_position('bottom')
      axxiell.tick_params(axis='x',reset=False,which='both',length=8,width=2)

    ## add labels.
    axwp.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=20)
    if(np.fabs(wppow) > 0.01):
      axwp.set_ylabel(r'$r_{\sigma}^{%.1f} w_p(r_{\sigma})$' % (wppow),fontsize=16)
    else:
      axwp.set_ylabel(r'$w_p(r_{\sigma}) \, [h^{-1} {\rm Mpc}]$',fontsize=20)

    axxiell.set_xlabel(r'$s \, [h^{-1} {\rm Mpc}]$',fontsize=20)
    if(np.fabs(xiellpow) > 0.01):
      axxiell.set_ylabel(r'$s^{%.1f} \xi_{\ell}(s)$' % (xiellpow),fontsize=20)
    else:
      axxiell.set_ylabel(r'$\xi_{\ell}(s)$',fontsize=20)

    axxi2d.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
    axxi2d.set_ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)


    return axxi2d, axwp, axxiell, fancydict

if __name__ == '__main__':
  print 'hi'

