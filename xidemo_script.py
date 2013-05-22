# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#### exercise in classes and matplotlib.

import xi2d
import wp
import xiell
import xi

reload(xi2d)
reload(xiell)
reload(wp)

# <codecell>

## correlation function measurements, computed elsewhere in a parallel C code.
fN = 'xifiles/data/collidedBR-collate-cmass-dr10v7-N-FBBRang_xigrid.butterfly'
fS = 'xifiles/data/collidedBR-collate-cmass-dr10v7-S-FBBRang_xigrid.butterfly'

## a theoretical prediction for xi in the absence of any peculiar velocities
fTreal = 'xifiles/theory/makeHJcatalogv2zspace2.cat.zspace-1.xiell.nbins900.butterfly'

xiN = xi2d.xi2d(fN)
xiS = xi2d.xi2d(fS)
xireal = xi2d.xi2d(fTreal)

# this function uses the filename associated with the object, and looks in the appropriate spot for the correct optimal weighting between the data sets.
# it also does sanity checks that the same binning was used, etc.
xiD = xiN+xiS

#overload print statement too give some relevant information about where the data came from, what the binning looks like, etc.
print xiD

# <codecell>

#see what raw correlation looks like.
ff = plt.figure(figsize=[6,6])
ii=plt.imshow(xiD.xi.reshape(xiD.n1d,xiD.n1d),extent=[xiD.rsig.min(),xiD.rsig.max(),xiD.rpi.max(),xiD.rpi.min()])
plt.xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
plt.ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
plt.colorbar(ii)
ff.savefig("crappybutterfly.png")

# <codecell>

# let's do a log mapping of xi so we can see more structure.

ff = plt.figure(figsize=[6,6])
rsigS, rpiS, xiS, n1dS = xiD.symmetrize()
xiremap = np.log10(xiS/xiD.xi.min())/np.log10(xiD.xi.max()/xiD.xi.min())
xiremap = xi2d.getflip(xiremap,n1dS)
xl=xiD.rsig.max()
yl=xiD.rpi.max()

#plt.imshow(xiremap.reshape(xiD.n1d,xiD.n1d),extent=[xiD.rsig.min(),xiD.rsig.max(),xiD.rpi.max(),xiD.rpi.min()])
ii=plt.imshow(xiremap.reshape(n1dS,n1dS),extent=[-xl,xl,-yl,yl])
plt.colorbar(ii)
plt.xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
plt.ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
ff.savefig("nicebutterfly.png")

# <codecell>

## a contour plot is easiest to work with, though.
print len(xiD.xi), xiD.xi[0]
ff, ax=xiD.makedensityplot()
xiD.addcontour(ax=ax)
ff.savefig("nicebutterflywcontours.png")

# <codecell>

## if there were no velocities distorting the line-of-sight (rpi) galaxy separations, contours would be round.

## make some custom contour levels between maximum of xi (at 0,0) and the border of the x-axis (~10).
myc = xiD.getclevels(ncontours=7,ximax=xiD.xiinterp(0,0)*0.9,ximin=xiD.xiinterp(9.5,0))

def plotDvsT(xdata,xth,cdata='k',cth='b'):
  fxi=plt.figure(figsize=[12,6])
  ax1=fxi.add_subplot(1,2,1)
  ax2=fxi.add_subplot(1,2,2)  
  xdata.makecontourplot(ax=ax1,color=cdata)
  xth.addcontour(ax=ax1,color=cth)
  xdata.makecontourplot(ax=ax2,span=[-10,10,-10,10],clevlist=myc,color=cdata)
  xth.addcontour(ax=ax2,color=cth,clevlist=myc)
  ax1.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
  ax1.set_ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
  ax2.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
  ax2.set_ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
  return fxi
    
myp=plotDvsT(xiD,xireal)

myp.savefig("datavsnovel.png")

# <codecell>

## this doesn't seem to run properly in ipython.
#%run animatef.py
import os
os.system('python animatef.py')
## need to install some package before I can save animations as movies.

# <codecell>

## before we compare to theory, we want to do some data compression steps.
## read in objects that keep track of xi (as shown above) along with various data compression variants.
## they also have methods to compute the chi2 for the various data compression variants.

## final(ish) data vector after removing small bias caused by "fiber collisions" and imperfect knowledge of the 
## redshift distribution of our sources.
xD=xi.xi(xidatfname='xifiles/data/dr10v7combinedFAKE.xidat')

## 3 different theoretical models with different behavior for the motions internal to their 
#host dark matter halos (i.e., the gravitational potential wells in which they are bound)

xT1=xi.xi(xidatfname='xifiles/theory/currbestfit_noihv.xidat')
xT2=xi.xi(xidatfname='xifiles/theory/currbestfit_nocenv.xidat')
xT3=xi.xi(xidatfname='xifiles/theory/currbestfit.xidat')

## the most common one (wp(rsig)) is to integrate along the line of sight (rpi)
tlist = [xT1.wp, xT2.wp, xT3.wp]
lbllist = ['no v', 'no cen v', 'best fit']
clist=['b','g','r']

fwp=plt.figure(figsize=[12,6])
ax1=fwp.add_subplot(1,2,1)
ax2=fwp.add_subplot(1,2,2)

wpD=xD.wp

lylist=[1,0] # logyopt
plist=[0,1]
axlist=[ax1,ax2]
fmtlist=[None,'ko']
slist=[[wpD.rsig.min()*0.9,wpD.rsig.max()*1.1,wpD.wp.min()*0.9,wpD.wp.max()*1.1],\
       [wpD.rsig.min()*0.9,wpD.rsig.max()*1.1,130,270]]

for (ax, p, s,ly,f) in zip(axlist, plist, slist,lylist,fmtlist):
  wpD.makeplot(color='k',ax=ax,rppow=p,span=s,lbl='data',logyopt=ly,fmt=f)
  for (cnt, t, lbl) in zip(range(len(tlist)), tlist, lbllist):
    t.addcurve(ax=ax,color=clist[cnt],rppow=p,lbl=lbl)
    if(p==0):
      print 'chi2 for %s: ' % (lbllist[cnt]),wpD.chi2(t)
    
ax1.legend()
ax1.plot([1,1],[0.1,2000],'k--')
ax2.plot([1,1],[0.1,2000],'k--')
ax1.text(0.2,80,'gravitationally bound',fontsize=16)
ax1.text(0.6,50,r'$\leftarrow$',fontsize=20)
ax1.text(1.1,10,r'$\rightarrow$ gravitational infall',fontsize=16)

fwp.savefig("wpcomparetheory.png")

# <codecell>

## money plot
#
## We came up with another data compression strategy that retains information on the velocity structure
## these are the monopole (spherical average) and quadrupole (anisotropy) moments of the correlation function
## when we include both velocity effects we get a good chi2!

wppow=1
xiellpow=1.25

myc = xD.xi2d.getclevels(ncontours=3,ximax=xD.xi2d.xiinterp(0,0)*0.9,ximin=xD.xi2d.xiinterp(9.5,0))
myc=[]

fig, a1, a2, a3, pset = xD.makefancyplot(sizescale=1,spanxiell=[0.2,35,-10,50],spanwp=[0.2,35,120,260],wppow=wppow,logyoptwp=0,xiellpow=xiellpow,color2='m',lbl='data')

#a1, a2, a3, pset = xD.makefancyplot(sizeratio=1.0,spanxiell=[0.2,35,-10,50],spanwp=[0.2,35,120,260],spanxi2d=[-10,10,-10,10],clevlist=myc,\
#wppow=wppow,logyoptwp=0,xiellpow=xiellpow,color2='g')

tlist = [xT1, xT2, xT3]
print 'number of degrees of freedom: ',len(xD.xiell.xilong)
for (cval, t, lbl) in zip(clist, tlist, lbllist):
  pset['color'] = pset['color2'] = cval
  pset['lbl'] = lbl
  t.addfancyplot(a1,a2,a3,pset)
  print 'chi2 for this theory is',xD.xiell.chi2(t.xiell)

a2.legend(loc=2)

fig.savefig("fancyfigv0.png")

# <codecell>

fxielldiff=plt.figure(figsize=[12,6])
ax1=fxielldiff.add_subplot(1,2,1)
ax2=fxielldiff.add_subplot(1,2,2)

elllist=[0,2]

axlist=[ax1,ax2]

for (ax, ell) in zip(axlist, elllist):
  for (cnt, t, lbl) in zip(range(len(tlist)), tlist, lbllist):
    #print 'yo',cnt,lbl
    #print t
    #print t.ndata
    xD.xiell.makediffplot(t.xiell,ell=ell,ax=ax,color=clist[cnt],lbl=lbl)
    if(ell==0):
      print 'chi2 for %s: ' % (lbllist[cnt]),xD.xiell.chi2(t.xiell)
    
#ax1.legend()
#ax2.legend()
ax1.axis([0.2,35,-5,5])
ax2.axis([0.2,35,-5,5])

fxielldiff.savefig("xielldiffcomparetheorytoerrors.png")

# <codecell>


