# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

try:
  reload(xismu)
except:
  import xismu

# <codecell>

mtaglist=['Mmin12.182Mmax12.483','Mmin12.484Mmax12.784','Mmin12.785Mmax13.085','Mmin13.086Mmax13.386','Mmin13.387Mmax17.000']
cnt=0
xilist=[]
myc=['k','b','g','r','c']

## this is how many mu bins we want; must be something that divides evenly into 200, at least for now.
getnmu = 4
rebinfac = 200/getnmu
for mtag in mtaglist:
  ffmt="/Users/bareid/work/montserratdata/forlilev3/A%02d_0.6452.halos.zspace%d.Np.nbins9000.bin1sims_"+mtag
  ## tmp input.
  xx=xismu.initavg(ffmt,zspacelist=[0,1,2],simlist=[0,2,3,4,5])
  ## final input.
  #xx=initavg(ffmt,[0,1,2],np.arange(0,20,1))
  ##rebin by a factor of rebinfac
  xilist.append(xx.rebin(1,rebinfac))
  outfname='avgmbin%ddown.smuxi' % (cnt+1)
  xilist[cnt].printsmuxi(outfname)
  if(cnt == 0):
    axlist=xilist[cnt].makexismuplot(panelopt=0,logxopt=0,color=myc[cnt])
  else:
    for jj in range(xilist[cnt].nmu):
      xilist[cnt].addcurve(axlist[0], jj,color=myc[cnt])

#   outplot = "mbin%d.smuxi.png" % cnt
#  xilist[cnt].fig.savefig(outplot)
  cnt += 1

plt.show()

# <codecell>

lbllist=[]
for muval in xilist[0].mu1d:
  ll = r'$\mu = %.3f$' % (muval)
  lbllist.append(ll)

print lbllist
## least massive bin.
aa=xilist[0].makexismuplot(panelopt=0,logxopt=0,color=['b','g','r','c'],lbl=lbllist)
plt.legend()
plt.title('low mass bin')
plt.show()

#most massive bin
aa=xilist[-1].makexismuplot(panelopt=0,logxopt=0,color=['b','g','r','c'],lbl=lbllist)
plt.legend()
plt.title('high mass bin')
plt.show()


# <codecell>


