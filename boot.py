import numpy as np
import sys
import fitsio
from scipy.integrate import quad

pixtype = [('PID','int'),('idec','int'),('ira','int'),('ramin','float'),('ramax','float'),('decmin','float'),('decmax','float'),('nR','int'),('nRtot','int'),('nD','int'),('PIDnew','int'),('NorS','int')]

def readsortradecwgtfromtxt(fname,targcat=0):
  if targcat == 1:
    ra, dec, wgt = np.loadtxt(fnameR,unpack=True)
  else:
    ra, dec, wgt = np.loadtxt(fnameR,unpack=True,usecols=[0,1,3])

  ra = ra - 90.0
  xx = np.where(ra < 0.0)[0]
  ra[xx] += 360.0
  del xx

  ii = dec.argsort()
  ra = ra[ii]
  dec = dec[ii]
  wgt = wgt[ii]
  return ra,dec,wgt


## copying ~/boss/bootstrapclean/bootall.py, taking out lots of hard coding but keep method the same. 
class bootpix:
  ## don't really need an init?
  def __init__(self,nsub):
    """
    Creates an empty bootpix.
    """
    self.nsub = -1
    self.pixlist = None

  def makeregions(self,nsub,maskareaNSdeg,decoffsetdeg,ddir='/home/howdiedoo/boss/mksamplecatslatestdr12/',sample='cmass',runtag='dr12v4',cattag='Reid-targ',sysopt=0,fkpopt=0,threshhold = 0.9):
  """
  nsub = number of subarrays we want to divide into.
  maskarea is read out of mask files; we want sum of N and S.
  decoffsetdeg is the value of dec to start the boot strap region constructions.  First elt is N, second is S.
  decoffsetdeg must be a numpy array!
  For dr10, we used 
  decoffsetsdeg = np.array([-3.0,-3.5]) 
  totarea = 5172.9866 + 1429.8623 = 6602.8489
  For DR12 -- 
  hmmm.
  """

  assert sysopt == 0 or sysopt == 1
  assert fkpopt == 0 or fkpopt == 1
  if sampletag == 'cmass':
    assert sysopt == 1
  else:
    assert sysopt == 0

  catappend = ''
  if sysopt == 1:
    catappend = catappend + '-wsys'
  if fkpopt == 1:
    catappend = catappend + '-wfkp'

  targcat = 0
  if re.search('targ',cattypetag):
    targcat = 1  # no redshifts!
  assert targcat == 1 # I think we want to stick to target catalogs for this?  Doesn't really matter much though.
  ## if we don't enforce this, then want to apply zcuts before computing the regions?

  degtorad = np.pi/180.
  NStag = ['N','S']

  ## record calling choices; print them in the pixel file?
  self.nsub = nsub
  self.marea = maskareaNS*degtorad**2
  assert len(decoffsetdeg) == 2
  self.decoffsets = decoffsetdeg*degtorad
  self.threshhold = threshhold
  self.fnameout = ddir + sampletag + '-' + runtag + '-' + cattypetag + catappend + boottag + '.Nsub%03d' % (self.nsub)

  ## set up default pixel size.
  self.ddec = (self.marea/self.nsub)**0.5
  nsidedec = int(np.floor(np.pi/ddec)+1)

  ## create structure to store pixels.
  pixmax = self.nsub*10
  pixlist = np.zeros(pixmax,dtype=pixtype)
  pixlistfinal = np.zeros(nsub,dtype=pixtype)
  mypix = 0 #location in pixlist.

  ## stores N and S pixel lists.
  pixRlist = []
  pixDlist = []
  ralist = []
  declist = []

  ## create the pixels in N and S independently.
  for NorS in [0,1]:
    first = 0

    myoffset = ((self.decoffsets[NorSopt]+0.5*np.pi)/ddec - np.floor((self.decoffsets[NorSopt]+0.5*np.pi)/ddec))*ddec
    if(myoffset > 0.5*ddec):
      assert np.fabs(myoffset) > np.fabs(myoffset-ddec)
      myoffset = myoffset-ddec

    NStag = NSlist[NorS]
    fnameD = ddir + sampletag + '-' + runtag + '-' + NStag + '-' + cattypetag + catappend + '.dat.txt'
    fnameR = ddir + sampletag + '-' + runtag + '-' + NStag + '-' + cattypetag + catappend + '.ran.txt'

    ra,dec,wgt = readsortradecwgtfromtxt(fnameR,targcat=targcat):
    raD,decD,wgtD = readsortradecwgtfromtxt(fnameD,targcat=targcat):

    ## save pixel index for every random/data object.
    pixR = np.zeros(len(ra),dtype='int')-1
    pixD = np.zeros(len(raD),dtype='int')-1

    ## loop over dec pixels.
    for i in range(nsidedec):
      decmin = -0.5*np.pi + i*ddec + myoffset
      decmax = -0.5*np.pi + (i+1)*ddec + myoffset
      cdec = np.cos(-0.5*np.pi + (i+0.5)*ddec + myoffset)
      ## do integral exactly.
      dint = quad(lambda x: np.cos(x), decmin, decmax)[0]
      dra = ddec**2/dint # size of pixel in ra direction.

      if(decmax <= dec[0] and decmax <= decD[0]):
        continue
      if(decmin > dec[-1] and decmin > decD[-1]):
        continue
      if(first == 0):
        ilow = 0
        ihigh = 0
        ilowD = 0
        ihighD = 0
        first = 1
      else:
        ilow = ihigh
        ilowD = ihighD

      print 'working on dec row',i,decmin,decmax

      while(dec[ihigh] < decmax and ihigh < len(dec)-1):
        ihigh += 1

      while(decD[ihighD] < decmax and ihighD < len(decD)-1):
        ihighD += 1

      if(ihigh == ilow):
        continue

      ## ihigh is not in the dec range, make loop go from ilow to ihigh-1.
      striperamin = min(ra[ilow:ihigh])
      striperamax = max(ra[ilow:ihigh])
      if(ihighD > ilowD):
        striperamin = min(min(ra[ilow:ihigh]),min(raD[ilowD:ihighD]))
        striperamax = max(max(ra[ilow:ihigh]),max(raD[ilowD:ihighD]))

      assert (dec[ilow:ihigh] < decmax).all()
      assert (dec[ilow:ihigh] >= decmin).all()
      assert (decD[ilowD:ihighD] < decmax).all()
      assert (decD[ilowD:ihighD] >= decmin).all()
      if(ilow > 0):  assert dec[ilow-1] < decmin
      if(ihigh < len(dec)-1):  assert dec[ihigh+1] >= decmax
      if(ilowD > 0):  assert decD[ilowD-1] < decmin
      if(ihighD < len(decD)-1):  assert decD[ihighD+1] >= decmax

      ## now optimize the centering of pixels in ra.
      ## screw it, let's do that later.
      nra = int(np.floor((striperamax-striperamin)/dra))+1
      extra = (nra*dra - (striperamax-striperamin))
      assert extra >= 0.
      rastart = striperamin - extra*0.5
      ## add an extra box, just because of the offseting.
      nra = nra + 1
      ## optimize on ra placement on the stripe -- want the most boxes "close to full"
      ## which we'll define as the most boxes above 90% full.  and if we need a tie breaker, the total randoms in those boxes should be highest.
      ## do this the dumb way; scan through a bunch of ra offsets.

      nobjlist = np.zeros(nra)
      ngoodbest = 0
      sumtotbest = 0
      bestraoff = -1000.
      for ioff in range(-50,51):
        myraoff = ioff/100.*dra
        for j in range(nra):
          ramin = rastart + j*dra + myraoff
          ramax = rastart + (j+1)*dra + myraoff
          if(j == nra-1):
            assert ramax > striperamax
          obj = np.where((ra[ilow:ihigh] >= ramin) & (ra[ilow:ihigh] <= ramax))[0]
          nobj = len(obj)
          nobjlist[j] = nobj
          ngood, sumtot = nobjtest(nobjlist,threshhold)
          if((ngood > ngoodbest) or (ngood == ngoodbest and sumtot > sumtotbest)):
            ngoodbest = ngood
            sumtotbest = sumtot
            bestraoff = myraoff

      assert bestraoff > -999

      ## go through once more with the bestraoff.
      for j in range(nra):
        ramin = rastart + j*dra + bestraoff
        ramax = rastart + (j+1)*dra + bestraoff
        if(j == nra-1):
          assert ramax > striperamax
        obj = np.where((ra[ilow:ihigh] >= ramin) & (ra[ilow:ihigh] <= ramax))[0]
        objD = np.where((raD[ilowD:ihighD] >= ramin) & (raD[ilowD:ihighD] <= ramax))[0]
        nobj = len(obj)
        nobjD = len(objD)

        if(nobj > 0):
          pixR[ilow:ihigh][obj] = mypix
        if(nobjD > 0):
          pixD[ilowD:ihighD][objD] = mypix
        if(nobj > 0 or nobjD > 0):
          ##ofpraw.write('%d %d %d %e %e %e %e %d %d\n' % (mypix,i, j, ramin, ramax, decmin, decmax, nobj, nobjD))
          pixlist['PID'][mypix] = mypix
          pixlist['idec'][mypix] = i
          pixlist['ira'][mypix] = j
          pixlist['ramin'][mypix] = ramin
          pixlist['ramax'][mypix] = ramax
          pixlist['decmin'][mypix] = decmin
          pixlist['decmax'][mypix] = decmax
          pixlist['nR'][mypix] = nobj
          pixlist['nRtot'][mypix] = nobj
          pixlist['nD'][mypix] = nobjD
          pixlist['PIDnew'][mypix] = -1
          pixlist['NorS'][mypix] = NorSopt
          mypix += 1
          assert mypix < pixmax

    #store pixR vals separately for N and S.
    pixRlist.append(pixR)
    pixDlist.append(pixD)
    ralist.append(ra)
    declist.append(dec)

  ## now we'll regroup these pixels to make a more uniform distribution of number of objects.
  npix = mypix # number of pixels we have.
  pixlist = pixlist[:npix]

  ## let's make a histogram of distributions in pixR, pixD
  ## skip for now!
  if 0==0:
    for i in range(len(pixlist)):
      xRN = np.where(pixRlist[0][:] == pixlist['PID'][i])[0]
      xDN = np.where(pixDlist[0][:] == pixlist['PID'][i])[0]
      xRS = np.where(pixRlist[1][:] == pixlist['PID'][i])[0]
      xDS = np.where(pixDlist[1][:] == pixlist['PID'][i])[0]
      print pixlist['nR'][i], len(xRN)+len(xRS), pixlist['nD'][i], len(xDN)+len(xDS)

  pixlist.sort(order='nRtot')
  print 'before concats, this randoms in top nsub:',(pixlist['nRtot'][-nsub:]).sum()/float((pixlist['nRtot'][:]).sum())
  mymed=0.5*(pixlist['nRtot'][-nsub/2]+pixlist['nRtot'][-nsub/2-1])
  print 'median occupancy',mymed
  print 'fractional var about median',(((pixlist['nRtot'][-nsub:]-mymed)**2).sum()/float(len(pixlist['nRtot'][-nsub:])))**0.5/mymed

  ## range we want groups of pixels to fall into.
  nRtarglow = pixlist['nRtot'][-50]
  diff = (pixlist['nRtot'][-1]-pixlist['nRtot'][-50])
  nRtarghigh = pixlist['nRtot'][-1] + diff
  print 'targs',nRtarglow,nRtarghigh,pixlist['nRtot'][-1]


  ## back to sorted on idec.
  pixlist.sort(order=('NorS','idec','ramin'))
  ## look for groups of 2 first -- only neighbors in dec.  later we can divide up small regions across different idecs more easily.

  for ii in range(len(pixlist)):
    #check for possibility of grouping.  if already grouped, continue
    if(pixlist['PIDnew'][ii] != -1):
      continue
    #find neighbors with same idec.
    ## real neighbors.
    nbest = -1
    ibest = -1
    for inbr in [ii-1,ii+1]:
      if inbr < 0 or inbr >= len(pixlist):
        continue
      if(pixlist['NorS'][ii] == pixlist['NorS'][inbr] and pixlist['idec'][ii] == pixlist['idec'][inbr]):
        ## one of RA borders should be shared.
        ## unless htere is a gap -- happens in N in current pixelization.
        ## don't join non-adjacent pixels.
        if(not(np.fabs(pixlist['ramin'][ii] - pixlist['ramax'][inbr]) < 2.0e-6 or np.fabs(pixlist['ramax'][ii] - pixlist['ramin'][inbr]) < 2.0e-6)):
          continue
        #nsum = pixlist['nR'][ii] + pixlist['nR'][inbr]
        ## switch to using nRtot, in case inbr has already joined with another.
        nsum = pixlist['nRtot'][ii] + pixlist['nRtot'][inbr]
        if(nsum >= nRtarglow and nsum <= nRtarghigh):
          nbest = nsum
          ibest = inbr
    if(ibest != -1): # found a match.
      if(pixlist['nRtot'][ii] > pixlist['nRtot'][ibest]):
        ibig = ii
        ismall = ibest
      else:
        ibig = ibest
        ismall = ii

      ## adding this so hopefully code will group threes correctly.
      if(pixlist['PIDnew'][ibest] != -1):
        ibig = ibest
        ismall = ii

      assert nbest == pixlist['nRtot'][ii] + pixlist['nRtot'][ibest]
      pixlist['nRtot'][ibig] = nbest
      pixlist['nRtot'][ismall] = -nbest
      pixlist['PIDnew'][ibig] = pixlist['PID'][ibig]
      pixlist['PIDnew'][ismall] = pixlist['PID'][ibig]

      ## assin new pixel value in the pixR and pixD arrays too
      xRN = np.where(pixRlist[0][:] == pixlist['PID'][ismall])[0]
      xDN = np.where(pixDlist[0][:] == pixlist['PID'][ismall])[0]
      xRS = np.where(pixRlist[1][:] == pixlist['PID'][ismall])[0]
      xDS = np.where(pixDlist[1][:] == pixlist['PID'][ismall])[0]
      if (len(xRN)+len(xRS)) != pixlist['nR'][ismall]:
        print 'wuuttR?',len(xRN),len(xRS),pixlist['nR'][ismall],ismall,pixlist['PID'][ismall]
        print 'maybe this pixel merged twice?',pixlist['PID'][ismall]
      ### this won't be true if merged two pixels ??
      #assert (len(xRN)+len(xRS)) == pixlist['nR'][ismall]
      if(len(xRN) > 0):
        pixRlist[0][xRN] = pixlist['PID'][ibig]
      if(len(xRS) > 0):
        pixRlist[1][xRS] = pixlist['PID'][ibig]

      if (len(xDN)+len(xDS)) != pixlist['nD'][ismall]:
        print 'wuuttD?',len(xDN),len(xDS),pixlist['nD'][ismall],ismall,pixlist['PID'][ismall]
        print 'maybe this pixel merged twice?',pixlist['PID'][ismall]
      #assert (len(xDN)+len(xDS)) == pixlist['nD'][ismall]
      if(len(xDN) > 0):
        pixDlist[0][xDN] = pixlist['PID'][ibig]
      if(len(xRS) > 0):
        pixDlist[1][xDS] = pixlist['PID'][ibig]

      print 'joined these two',pixlist['nR'][ibig],pixlist['nR'][ismall],nsum, pixlist['PID'][ismall], pixlist['PID'][ibig]
      print pixlist[ibig]
      print pixlist[ismall]

  # now go through pixels that are unassigned and tack them to nearest dec neighbor, draw new ra lines if necessary.
  ## need to code.

  ## now get top nsub
  pixlist.sort(order='nRtot')
  print 'after concats, this randoms in top nsub:',(pixlist['nRtot'][-nsub:]).sum()/float((pixlist['nR'][:]).sum())
  mymed=0.5*(pixlist['nRtot'][-nsub/2]+pixlist['nRtot'][-nsub/2-1])
  print 'median occupancy',mymed
  print 'fractional var about median',(((pixlist['nRtot'][-nsub:]-mymed)**2).sum()/float(len(pixlist['nRtot'][-nsub:])))**0.5/mymed

  ## optionally print out the raw pixlist(?)

  ## now fill in pixlistfinal from the original pixlist.
  for ii in range(-nsub,0):
    if(pixlist['nR'][ii] > pixlist['nRtot'][ii]):
      sys.exit(1)
    if(pixlist['nR'][ii] < pixlist['nRtot'][ii]):
      xxx = np.where(pixlist['PIDnew'] == pixlist['PIDnew'][ii])[0]
      assert ((pixlist['PIDnew'][xxx] == pixlist['PIDnew'][ii]) | (xxx < len(pixlist)-nsub)).all()
      assert (pixlist['decmin'][xxx] == pixlist['decmin'][ii]).all()
      assert (pixlist['decmax'][xxx] == pixlist['decmax'][ii]).all()
      myramin = (pixlist['ramin'][xxx]).min()
      myramax = (pixlist['ramax'][xxx]).max()
      mydecmin = (pixlist['decmin'][xxx]).min()
      mydecmax = (pixlist['decmax'][xxx]).max()

      pindx = nsub + ii
      
      pixlistfinal['PID'][pindx] = pindx
      pixlistfinal['idec'][pindx] = pixlist['idec'][ii]
      pixlistfinal['ira'][pindx] = pixlist['ira'][ii]
      pixlistfinal['ramin'][pindx] = myramin
      pixlistfinal['ramax'][pindx] = myramax
      ## ahh! this is what DR10 version did, but we want mydecmin/mydecmax!
      #pixlistfinal['decmin'][pindx] = pixlist['decmin'][ii]
      #pixlistfinal['decmax'][pindx] = pixlist['decmax'][ii]
      pixlistfinal['decmin'][pindx] = mydecmin
      pixlistfinal['decmax'][pindx] = mydecmax

      xRN = np.where(pixRlist[0][:] == pixlist['PIDnew'][ii])[0]
      xRS = np.where(pixRlist[1][:] == pixlist['PIDnew'][ii])[0]
      xDN = np.where(pixDlist[0][:] == pixlist['PIDnew'][ii])[0]
      xDS = np.where(pixDlist[1][:] == pixlist['PIDnew'][ii])[0]
      assert len(xRN) > 0 or len(xRS) > 0
      pixlistfinal['nR'][pindx] = len(xRN) + len(xRS)
      pixlistfinal['nRtot'][pindx] = len(xRN) + len(xRS)
      pixlistfinal['nD'][pindx] = len(xDN) + len(xDS)
      pixlistfinal['PIDnew'][pindx] = pindx
      pixlistfinal['NorS'][pindx] = pixlist['NorS'][ii]

    else:
      pixlistfinal[pindx] = pixlist[ii]
      xRN = np.where(pixRlist[0][:] == pixlist['PID'][ii])[0]
      xRS = np.where(pixRlist[1][:] == pixlist['PID'][ii])[0]
      assert len(xRN) > 0 or len(xRS) > 0

    if not (len(xRN) == pixlist['nRtot'][ii] or len(xRS) == pixlist['nRtot'][ii]):
      print 'this should be assert error!',len(xRN),len(xRS),pixlist['nRtot'][ii]
    assert len(xRN) == pixlist['nRtot'][ii] or len(xRS) == pixlist['nRtot'][ii]

  self.pixlist = pixlistfinal
  ## finished makeregions.

#pixtype = [('PID','int'),('idec','int'),('ira','int'),('ramin','float'),('ramax','float'),('decmin','float'),('decmax','float'),('nR','int'),('nRtot','int'),('nD','int'),('PIDnew','int'),('NorS','int')]
  def writeregions(self):
    ofp = open(self.fnameout,'w')
    ## write out input params.
    ofp.write('nsub: %d\n' % (self.nsub))
    ofp.write('marea: %e\n' % (self.marea))
    ofp.write('decoffsets: %e, %e\n' % (self.decoffsets[0],self.decoffsets[1]))
    ofp.write('threshhold: %e\n' % (self.threshhold))
    ofp.write('ddec: %e\n' % (self.ddec))
    p = self.pixlist #shorthand.
    for pi in range(self.nsub):
      ofp.write('%d %d %d %d %e %e %e %e\n' % (p['PID'][pi], p['NorS'][pi], p['idec'][pi], p['ira'][pi], p['ramin'][pi], p['ramax'][pi], p['decmin'][pi], p['decmax'][pi],p['nRtot'][pi]))
    ofp.close()
  
  def readregions(self,fname):
    aa = np.loadtxt(fname)

    self.nsub = len(aa[:,0])
    assert len(aa[0,:]) == 8
    self.pixlist = np.zeros(nsub,dtype=pixtype)
    self.fnameout = fname
    pdict = {'PID':0,'NorS':1,'idec':2,'ira':3,'ramin':4,'ramax':5,'decmin':6,'decmax':7,'nRtot':8}
    for k,v in pdict.iteritems():
      self.pixlist[k] = aa[:,v]
    self.pixlist['nR'][:] = -1
    self.pixlist['nD'][:] = -1
    self.pixlist['PIDnew'][:] = -1

    ifp.close()

  def plotregions(self,ax=None):
    """
    Plot the borders of the regions.
    Need to code!
    """
    pass

  def writesubcats(self):
    """
    Need to copy bootstrapclean/writesubcats200.py, pay attention to orphan option.
    """
    pass

  def plotsubcats(self,ax=None):
    """
    Plot subregions by different colors to visualize them.
    """
    pass

class bootsec:  ## also divide survey into bootstrap regions generated by joining sectors together up to some area.
  __init__(self):
    pass


