import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
import types 
import os
import sys
import xi

flist = np.arange(0.6,1.4001,0.05)
fpath1 = "../../work/montserratdata/boss/bethalexie/"
fpath2 = "../../work/montserratdata/boss/bethalexie/fits2data/"
foutpath = "xifiles/theory/"

## run this once to move data over and generate xidat files.
if(0==1):

  for fval in flist:
    if fval < 1.:
      #tag = '0p'+str(int(fval*1000))
      tag = '0p'+'%03d' % (int(fval*1000 + 0.5))
    else:
      #tag = '1p'+str(int((fval-1.)*1000+0.5))
      tag = '1p'+'%03d' % (int((fval-1.)*1000 + 0.5))
  
    f1 = "makeHJcatalogzspace2v2_varyf"+tag+".cat.zspace2.xiell.nbins900.butterfly"
    f2 = "makeHJcatalogzspace2v2_varyf"+tag+".xiell.0"
    f3 = "makeHJcatalogzspace2v2_varyf"+tag+".wp.0"

    mystr = 'cp %s%s %s' % (fpath1,f1,foutpath)
    os.system(mystr)

    mystr = 'cp %s%s %s' % (fpath2,f2,foutpath)
    os.system(mystr)

    mystr = 'cp %s%s %s' % (fpath2,f3,foutpath)
    os.system(mystr)
  
    outf = foutpath + "makeHJcatalogv2zspace2_varyf"+tag+".xidat"
    ofp = open(outf,'w')
    ofp.write('xi2d: %s%s\n' % (foutpath,f1))
    ofp.write('wp: %s%s\n' % (foutpath,f3))
    ofp.write('xiell: %s%s\n' % (foutpath,f2))
    ofp.close()
  sys.exit(1)

ff=plt.figure(figsize=[6,6])
ax = ff.add_subplot(1,1,1)
ax.set_xlabel(r'$r_{\sigma} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
ax.set_ylabel(r'$r_{\pi} \, [h^{-1} {\rm Mpc}]$',fontsize=16)
ax.set_title('infall velocities scaled by 0.6-1.4')

ims = []
for (cnt, fval) in zip(range(len(flist)), flist):
  if fval < 1.:
    #tag = '0p'+str(int(fval*1000))
    tag = '0p'+'%03d' % (int(fval*1000 + 0.5))
  else:
    #tag = '1p'+str(int((fval-1.)*1000+0.5))
    tag = '1p'+'%03d' % (int((fval-1.)*1000 + 0.5))

  outf = foutpath + "makeHJcatalogv2zspace2_varyf"+tag+".xidat"
  xif = xi.xi(xidatfname=outf)
#  print xif.wp
#  print xif.wp.rsig

  if(cnt == 0):
    clev=xif.xi2d.getclevels(ximin=xif.xi2d.xiinterp(20,0),ximax=xif.xi2d.xiinterp(0.,0.)*0.5,cratio=0.7)
    f_template = 'f = %.2f'
    f_text = ax.text(-25,25,'',transform=ax.transAxes)
#  print clev
  im = xif.xi2d.addcontour(ax=ax,clevlist=clev)
  if(np.fabs(fval - 1.) < 1.0e-3):
    xif.xi2d.addcontour(ax=ax,clevlist=clev,color='b')

  def setvisible(self,vis):
    for c in self.collections: c.set_visible(vis)
  im.set_visible = types.MethodType(setvisible,im,None)
  ### this doesn't seem to work, oh well!
  #f_text.set_text(f_template%fval) 
  im.axes = plt.gca() 
  ims.append([im]) 

ani = animation.ArtistAnimation(ff, ims, interval=400, blit=False, 
    repeat_delay=300) 

  #axxi2d, axwp, axxiell, pset = xif.makefancyplot()
plt.show()
