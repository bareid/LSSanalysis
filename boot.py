import numpy as np
import matplotlib.pyplot as plt
import sys
import fitsio
from scipy.integrate import quad
import time
import re

def nobjtest(nobjlist,threshhold):
  xx = np.where(nobjlist > threshhold*nobjlist.max())[0]
  ngood = 0
  sumtot = 0
  if(len(xx) > 0):
    ngood = len(xx)
    sumtot = (nobjlist[xx]).sum()
  return ngood, sumtot

pixtype = [('PID','int'),('idec','int'),('ira','int'),('ramin','float'),('ramax','float'),('decmin','float'),('decmax','float'),('nR','int'),('nRtot','int'),('nD','int'),('PIDnew','int'),('NorS','int')]
pixtypeout = [('PID','int'),('NorS','int'),('idec','int'),('ira','int'),('ramin','float'),('ramax','float'),('decmin','float'),('decmax','float'),('nRtot','int')]

def readsortradecwgtfromfits(fitsfname,DorR,decsortopt,cpopt,sysopt,fkpopt=0,ramin=-1000.,ramax=1000.,decmin=-1000.,decmax=1000.):

  dd = fitsio.read(fitsfname)
  t1 = time.time()
  xx = np.where((dd['RA'] >= ramin) & (dd['RA'] <= ramax) & (dd['DEC'] >= decmin) & (dd['DEC'] <= decmax))[0]
  t2 = time.time()
  print 'timing where statement:',t2-t1,'seconds'
  print 'keeping',len(xx),'out of',len(dd)
  dd = dd[xx]

  ra = np.zeros(len(dd))
  dec = np.zeros(len(dd))
  z = np.zeros(len(dd))
  ra[:] = dd['RA'][:]
  dec[:] = dd['DEC'][:]
  try:
    z[:] = dd['Z'][:]
  except:
    z[:] = -1.0

  wgtvec = np.zeros(len(dd))

  if cpopt == 0 or DorR == 1:
    wgtvec = np.zeros(len(dd)) + 1.
  if cpopt == 1 and DorR == 0:
    wgtvec = dd['WEIGHT_NOZ']
  if cpopt == 2 and DorR == 0:
    wgtvec = (dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0)
  if sysopt == 1 and DorR == 0:
    wgtvec = wgtvec * dd['WEIGHT_SYSTOT']
  if fkpopt == 1:
    wgtvec = wgtvec * dd['WEIGHT_FKP']

  del dd
  ra = ra - 90.0
  xx = np.where(ra < 0.0)[0]
  ra[xx] += 360.0
  del xx

  if decsortopt != 0:
    ii = dec.argsort()
    ra = ra[ii]
    dec = dec[ii]
    z = z[ii]
    wgtvec = wgtvec[ii]

  ra = ra*np.pi/180.
  dec = dec*np.pi/180.
  return ra,dec,z,wgtvec


## copying ~/boss/bootstrapclean/bootall.py, taking out lots of hard coding but keep method the same. 
class bootpix:
  ## don't really need an init?
  def __init__(self):
    """
    Creates an empty bootpix.
    """
    self.nsub = -1
    self.pixlist = None
    self.setup = 0 # is look-up table set up?  Only set it up if necessary.

  def makeregions(self,nsub,maskareaNSdeg,decoffsetdeg,boottag, NorSopt=2, ddir='/home/howdiedoo/boss/mksamplecatslatestdr12/',sampletag='cmass',runtag='dr12v4',cattypetag='Reid-targ',sysopt=1,fkpopt=0,threshhold = 0.9,raminIN=-1000.,ramaxIN=1000.,decminIN=-1000.,decmaxIN=1000.):
    """
    nsub = number of subarrays we want to divide into.
    maskarea is read out of mask files; we want sum of N and S.
    decoffsetdeg is the value of dec to start the boot strap region constructions.  First elt is N, second is S.
    decoffsetdeg must be a numpy array!
    For dr10, we used 
    decoffsetsdeg = np.array([-3.0,-3.5]) 
    totarea = 5172.9866 + 1429.8623 = 6602.8489
    For DR12 -- hmmm.
    For testing, can just cut out a small patch using ramin/max, decmin/max.
    """
    assert NorSopt >= 0 and NorSopt <= 2
    if NorSopt == 2:
      NorSlist = [0,1]
    else:
      NorSlist = [NorSopt]  
  
    
  
    ## default, for regular Reid catalogs.
    cpopt = 2
  
    assert sysopt == 0 or sysopt == 1
    assert fkpopt == 0 or fkpopt == 1
    if sampletag == 'cmass':
      assert sysopt == 1
    else:
      assert sysopt == 0
   
    targcat = 0
    if re.search('targ',cattypetag):
      targcat = 1  # no redshifts!
      cpopt = 0
    if re.search('ang',cattypetag):
      cpopt = 1
  
    assert targcat == 1 # I think we want to stick to target catalogs for this?  Doesn't really matter much though.
    ## if we don't enforce this, then want to apply zcuts before computing the regions?
  
    degtorad = np.pi/180.
    NSlist = ['N','S']
  
    ## record calling choices; print them in the pixel file?
    self.nsub = nsub
    self.marea = maskareaNSdeg*degtorad**2
    assert len(decoffsetdeg) == 2
    self.decoffsets = decoffsetdeg*degtorad
    self.threshhold = threshhold
    self.fnameout = ddir + sampletag + '-' + runtag + '-' + cattypetag + '-' + boottag + '.Nsub%04d' % (self.nsub)
  
    ## set up default pixel size.
    self.ddec = (self.marea/self.nsub)**0.5
    nsidedec = int(np.floor(np.pi/self.ddec)+1)

    print 'pixel size is ',self.ddec,self.ddec/degtorad,(self.ddec)**2/degtorad**2,maskareaNSdeg/(self.ddec)**2*degtorad**2
  
    ## create structure to store pixels.
    pixmax = self.nsub*1000
    pixlist = np.zeros(pixmax,dtype=pixtype)
    pixlistfinal = np.zeros(nsub,dtype=pixtype)
    mypix = 0 #location in pixlist.
  
    ## stores N and S pixel lists.
    pixRlist = []
    pixDlist = []
    ralist = []
    declist = []
  
    ## create the pixels in N and S independently.
    for NorS in NorSlist:
      first = 0
  
      myoffset = ((self.decoffsets[NorS]+0.5*np.pi)/self.ddec - np.floor((self.decoffsets[NorS]+0.5*np.pi)/self.ddec))*self.ddec
      if(myoffset > 0.5*self.ddec):
        assert np.fabs(myoffset) > np.fabs(myoffset-self.ddec)
        myoffset = myoffset-self.ddec
  
      NStag = NSlist[NorS]
      fnameD = ddir + sampletag + '-' + runtag + '-' + NStag + '-' + cattypetag + '.dat.fits'
      fnameR = ddir + sampletag + '-' + runtag + '-' + NStag + '-' + cattypetag + '.ran.fits'
  
      t1 = time.time()
      raD,decD,zD,wgtD = readsortradecwgtfromfits(fnameD,DorR=0,decsortopt=1,cpopt=cpopt,sysopt=sysopt,ramin = raminIN, ramax=ramaxIN, decmin = decminIN, decmax = decmaxIN)
      t2 = time.time()
      print 'read for data took',t2-t1,'seconds'
      ra,dec,z,wgt = readsortradecwgtfromfits(fnameR,DorR=1,decsortopt=1,cpopt=cpopt,sysopt=sysopt,ramin = raminIN, ramax=ramaxIN, decmin = decminIN, decmax = decmaxIN)
      t3 = time.time()
      print 'read for data took',t3-t1,'seconds'
      ## save pixel index for every random/data object.
      pixR = np.zeros(len(ra),dtype='int')-1
      pixD = np.zeros(len(raD),dtype='int')-1

  
      ## loop over dec pixels.
      for i in range(nsidedec):
        decmin = -0.5*np.pi + i*self.ddec + myoffset
        decmax = -0.5*np.pi + (i+1)*self.ddec + myoffset
        cdec = np.cos(-0.5*np.pi + (i+0.5)*self.ddec + myoffset)
        ## do integral exactly.
        dint = quad(lambda x: np.cos(x), decmin, decmax)[0]
        dra = self.ddec**2/dint # size of pixel in ra direction.
  
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
            pixlist['NorS'][mypix] = NorS
            mypix += 1
            assert mypix < pixmax
  
      #store pixR vals separately for N and S.
      pixRlist.append(pixR)
      pixDlist.append(pixD)
      ralist.append(ra)
      declist.append(dec)

    if mypix < nsub:
      print 'did not generate nsub subregions.  Adjust area?'
      return -1

    ## now we'll regroup these pixels to make a more uniform distribution of number of objects.
    npix = mypix # number of pixels we have.
    pixlist = pixlist[:npix]
  
    ## let's make a histogram of distributions in pixR, pixD
    ## skip for now! doesn't work right now because pixlists might only have one element.
    if 0==0:
      for i in range(len(pixlist)):
        xRN = np.where(pixRlist[0][:] == pixlist['PID'][i])[0]
        xDN = np.where(pixDlist[0][:] == pixlist['PID'][i])[0]
        if len(pixRlist) > 1:
          xRS = np.where(pixRlist[1][:] == pixlist['PID'][i])[0]
          xDS = np.where(pixDlist[1][:] == pixlist['PID'][i])[0]
        else:
          xRS = np.array([])
          xDS = np.array([])

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
        if len(pixRlist) > 1:
          xRS = np.where(pixRlist[1][:] == pixlist['PID'][ismall])[0]
          xDS = np.where(pixDlist[1][:] == pixlist['PID'][ismall])[0]
        else:
          xRS = np.array([])
          xDS = np.array([])

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
        xDN = np.where(pixDlist[0][:] == pixlist['PIDnew'][ii])[0]
        if len(pixRlist) > 1:
          xRS = np.where(pixRlist[1][:] == pixlist['PIDnew'][ii])[0]
          xDS = np.where(pixDlist[1][:] == pixlist['PIDnew'][ii])[0]
        else:
          xRS = np.array([])
          xDS = np.array([])

        assert len(xRN) > 0 or len(xRS) > 0
        pixlistfinal['nR'][pindx] = len(xRN) + len(xRS)
        pixlistfinal['nRtot'][pindx] = len(xRN) + len(xRS)
        pixlistfinal['nD'][pindx] = len(xDN) + len(xDS)
        pixlistfinal['PIDnew'][pindx] = pindx
        pixlistfinal['NorS'][pindx] = pixlist['NorS'][ii]
  
      else:
        pindx = nsub + ii
        pixlistfinal[pindx] = pixlist[ii]
        pixlistfinal['PID'][pindx] = pindx
        xRN = np.where(pixRlist[0][:] == pixlist['PID'][ii])[0]
        if len(pixRlist) > 1:
          xRS = np.where(pixRlist[1][:] == pixlist['PID'][ii])[0]
        else:
          xRS = np.array([])

        assert len(xRN) > 0 or len(xRS) > 0
  
      if not (len(xRN) == pixlist['nRtot'][ii] or len(xRS) == pixlist['nRtot'][ii]):
        print 'this should be assert error!',len(xRN),len(xRS),pixlist['nRtot'][ii]
      assert len(xRN) == pixlist['nRtot'][ii] or len(xRS) == pixlist['nRtot'][ii]
  
    self.pixlist = pixlistfinal
    ## finished makeregions.

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
      ofp.write('%d %d %d %d %e %e %e %e %d\n' % (p['PID'][pi], p['NorS'][pi], p['idec'][pi], p['ira'][pi], p['ramin'][pi], p['ramax'][pi], p['decmin'][pi], p['decmax'][pi],p['nRtot'][pi]))
    ofp.close()
  

  def readregions(self,fname):
    ifp = open(fname,'r')
    line = ifp.readline()
    ## read in parameters about the subregions.
    if not re.search('^nsub:',line):
      print 'missing nsub line'
      return None
    self.nsub = int(line.split(':')[1])
    
    line = ifp.readline()
    if not re.search('^marea:',line):
      print 'missing marea line'
      return None
    self.marea = float(line.split(':')[1])

    line = ifp.readline()
    if not re.search('^decoffsets:',line):
      print 'missing decoffsets line'
      return None
    self.decoffsets = np.array([float(line.split(':')[1].split(',')[0]), float(line.split(':')[1].split(',')[1])])

    line = ifp.readline()
    if not re.search('^threshhold:',line):
      print 'missing threshhold line'
      return None
    self.threshhold = float(line.split(':')[1])

    line = ifp.readline()
    if not re.search('^ddec:',line):
      print 'missing ddec line'
      return None
    self.ddec = float(line.split(':')[1])
    self.fnameout = fname
    #print self.nsub, self.marea, self.decoffsets, self.threshhold, self.ddec, self.fnameout
    ifp.close()
    aa = np.loadtxt(fname,skiprows=5,dtype=pixtypeout)

    assert(self.nsub == len(aa))
    self.pixlist = np.zeros(self.nsub,dtype=pixtype)
    pdict = {'PID':0,'NorS':1,'idec':2,'ira':3,'ramin':4,'ramax':5,'decmin':6,'decmax':7,'nRtot':8}
    for k,v in pdict.iteritems():
      self.pixlist[k] = aa[k]
    self.pixlist['nR'][:] = -1
    self.pixlist['nD'][:] = -1
    self.pixlist['PIDnew'][:] = -1
    self.setup = 0 # is look-up table set up?  Only set it up if necessary.

  def plotregions(self,ax=None,color='k'):
    """
    Plot the borders of the regions.
    Note ra has been shifted internally by 90 degrees throughout!
    We'll plot the shifted ra for simplicity.
    """
    degtorad = np.pi/180.
    if ax is None:
      ff = plt.figure(figsize=[6,6])
      ax=ff.add_subplot(1,1,1)
    else:
      ff = None

    offsetra = np.zeros(len(self.pixlist))
    ## nevermind, let's leave it offset coordinates.
    #offsetra = np.zeros(len(self.pixlist)) + np.pi/2.
    #xx = np.where(self.pixlist['ramin'] + offsetra > 2.*np.pi)[0]
    #offsetra[xx] -= 2.*np.pi

    for ii in range(self.nsub):
      r1 = (self.pixlist['ramin'][ii] + offsetra[ii])/degtorad
      r2 = (self.pixlist['ramax'][ii] + offsetra[ii])/degtorad
      d1 = self.pixlist['decmin'][ii]/degtorad
      d2 = self.pixlist['decmax'][ii]/degtorad

      ax.plot([r1,r1],[d1,d2],color=color)
      ax.plot([r2,r2],[d1,d2],color=color)
      ax.plot([r1,r2],[d1,d1],color=color)
      ax.plot([r1,r2],[d2,d2],color=color)

    return ff,ax

  #ra = ra - 90.0
  #xx = np.where(ra < 0.0)[0]
  #ra[xx] += 360.0

  def setuplookup(self):
    """
    Set up look-up table information after a pixelization is created.
    
    """
    ## additional parameters needed for look-up table.
    Np = np.where(self.pixlist['NorS'] == 0)[0]
    self.decminN = (self.pixlist['decmin'][Np]).min()
    xx = np.where(self.pixlist['decmin'][Np] == self.decminN)[0]
    self.idecminN = self.pixlist['idec'][Np[xx[0]]]
    ddecN = self.pixlist['decmax'][Np[0]] - self.pixlist['decmin'][Np[0]]
    assert (np.fabs(self.pixlist['decmax'][Np] - self.pixlist['decmin'][Np] - ddecN) < 2.0e-6).all()
    ## check sanity of idec formula.
  
    Sp = np.where(self.pixlist['NorS'] == 1)[0]
    self.decminS = (self.pixlist['decmin'][Sp]).min()
    xx = np.where(self.pixlist['decmin'][Sp] == self.decminS)[0]
    self.idecminS = self.pixlist['idec'][Sp[xx[0]]]
    ddecS = self.pixlist['decmax'][Sp[0]] - self.pixlist['decmin'][Sp[0]]
    assert (np.fabs(self.pixlist['decmax'][Sp] - self.pixlist['decmin'][Sp] - ddecS) < 2.0e-6).all()
    assert np.fabs(ddecN - ddecS) < 2.0e-6 #these should be the same
    ## should already be assigned!
    #self.ddec = ddecN
    assert np.fabs(ddecN - self.ddec) < 2.0e-6
    assert np.fabs(ddecS - self.ddec) < 2.0e-6
    idecchk = np.array(np.floor((self.pixlist['decmin'][Sp] - self.decminS)/ddecS+0.5)+self.idecminS,dtype='int')
    assert (idecchk == self.pixlist['idec'][Sp]).all()
    idecchk = np.array(np.floor((self.pixlist['decmin'][Np] - self.decminN)/ddecN+0.5)+self.idecminN,dtype='int')
    assert (idecchk == self.pixlist['idec'][Np]).all()
    print 'passed idec checks!'

    ## create a dictionary with a list of all the pixels with (NorS, idec) keys.
    ## each entry contains a global (over idec) gramin, gramax, decmin, decmax, ramin[], ramax[], PID[] which contain ramin, ramax, PID for all the pixels in that idec.
    pixdict = {}
    for p in self.pixlist:
      key = (p['NorS'], p['idec'])
      if key in pixdict:
        for k in ['decmin','decmax']:
          assert p[k] == pixdict[key][k]
        pixdict[key]['gramin'] = min(pixdict[key]['gramin'], p['ramin'])
        pixdict[key]['gramax'] = max(pixdict[key]['gramax'], p['ramax'])
        pixdict[key]['ramin'].append(p['ramin'])
        pixdict[key]['ramax'].append(p['ramax'])
        pixdict[key]['PID'].append(p['PID'])
      else:
        # entry is itself a dictionary!
        entry = {}
        entry['decmin'] = p['decmin']
        entry['decmax'] = p['decmax']
        entry['gramin'] = p['ramin']
        entry['gramax'] = p['ramax']
        entry['ramin'] = [p['ramin']]
        entry['ramax'] = [p['ramax']]
        entry['PID'] = [p['PID']]
        ## compute dra in this region from ddec.
        dint = quad(lambda x: np.cos(x), p['decmin'], p['decmax'])[0]
        dra = self.ddec**2/dint # size of pixel in ra direction.
        entry['dra'] = dra
        pixdict[key] = entry

    self.pixdict = pixdict 
    self.setup = 1 #look up table set up now!

  def getidec(self,NorS,dec):
    """
    Return integer dec pixel index for value dec.
    Currently works on a single dec value at a time.
    Easily upgradeable to numpy arrays if I want.
    """
    if(NorS == 0):
      return int(np.floor((dec - self.decminN)/self.ddec)+self.idecminN)
    else:
      return int(np.floor((dec - self.decminS)/self.ddec)+self.idecminS)

  def getpix(self,NorS,ra,dec,orphanopt):
    """
    Return pixel value ra,dec falls into (or nearest pixel or orphanopt = 1) 
    Currently works on a single ra/dec pair of values at a time.
    Easily upgradeable to numpy arrays if I want.
    Note input ra/dec are in radians, with 90 degree offset in ra already applied!
    orphanopt = 0 only returns pixel values if ra/dec inside the pixel.
    orphanopt will return nearby pixels.
    """
    printon = 0

    idec = self.getidec(NorS,dec)
    key = (NorS,idec)
    myval = self.pixdict.get(key)

    radiffB = 10000.
    decdiffB = 10000.
    dB = 10000.
    pB = -1

    if myval is not None:
      if(not(dec >= myval['decmin'] and dec < myval['decmax'])):
        print 'this assert went off',dec,myval['decmin'],myval['decmax']
        #assert dec >= myval['decmin'] and dec < myval['decmax']  ## I had commented this line out, why?
      for i in range(len(myval['ramin'])):
        ## definitely inside a pixel
        if(ra >= myval['ramin'][i] and ra < myval['ramax'][i]):
          return myval['PID'][i]
        ## this guy is inside dec boundaries and within dra/2 of ra border, but not in a pixel.
        radiff = 100000.
        if(ra >= myval['ramin'][i] - myval['dra']*0.5 and ra < myval['ramin'][i]):
          radiff = (myval['ramin'][i] - ra)*np.cos(dec)
          assert radiff > 0.
          assert radiff <= myval['dra']*np.cos(dec)*0.5
  
        if(ra < myval['ramax'][i] + myval['dra']*0.5 and ra >= myval['ramax'][i]):
          radiff = (ra - myval['ramax'][i])*np.cos(dec)
          assert radiff > 0.
          assert radiff <= myval['dra']*np.cos(dec)*0.5
        if(radiff < radiffB):
          if(printon == 1):
            print 'yo beth matching',ra,dec,(myval['ramin'][i]-ra)/myval['dra']/0.5,(ra-myval['ramax'][i])/myval['dra']/0.5,myval['PID'][i],myval
          pB = myval['PID'][i]
          radiffB = radiff
          dB = radiff
          decdiffB = 0.

    if orphanopt == 0:
      return -1 ## don't group with nearby pixel if orphanopt == 0

    ## if we got here, the object doesn't fall into a pixel.  Is it near one in the dec direction?
    ilist = [idec-1,idec+1]
    for idnew in ilist:
      key = (NorS,idnew)
      myval = self.pixdict.get(key)
      if myval is not None:
        decdist = min(np.fabs(dec-myval['decmin']),np.fabs(dec-myval['decmax']))
        if(decdist/self.ddec > 0.5001): continue  ##001 just to make sure a dec halfway between decmin and decmax is fine.
        if(decdist >= dB): continue # closer to the pixel above.
        for i in range(len(myval['ramin'])):
          ## inside ra pixel, dec is closer than ra.
          if(ra >= myval['ramin'][i] and ra < myval['ramax'][i]):
            pB = myval['PID'][i]
            decdiffB = decdist
            radiffB = 0.
            dB = decdist
          ## check near corners.
          radiff = 100000.
          if(ra >= myval['ramin'][i] - myval['dra']*0.5 and ra < myval['ramin'][i]):
            radiff = (myval['ramin'][i] - ra)*np.cos(dec)
            assert radiff > 0.
  
          if(ra < myval['ramax'][i] + myval['dra']*0.5 and ra >= myval['ramax'][i]):
            radiff = (ra - myval['ramax'][i])*np.cos(dec)
            assert radiff > 0.
          dtot = (radiff**2 + decdist**2)**0.5
          if(dtot < dB):
            pB = myval['PID'][i]
            radiffB = radiff
            decdiffB = decdist
            dB = dtot
  
    # if pB has been assigned, return that.  otherwise return -1
    return pB


#pixtype = [('PID','int'),('idec','int'),('ira','int'),('ramin','float'),('ramax','float'),('decmin','float'),('decmax','float'),('nR','int'),('nRtot','int'),('nD','int'),('PIDnew','int'),('NorS','int')]
  #def makeregions(self,nsub,maskareaNSdeg,decoffsetdeg,boottag, NorSopt=2, ddir='/home/howdiedoo/boss/mksamplecatslatestdr12/',sampletag='cmass',runtag='dr12v4',cattypetag='Reid-targ',sysopt=1,fkpopt=0,threshhold = 0.9,raminIN=-1000.,ramaxIN=1000.,decminIN=-1000.,decmaxIN=1000.):

  def writesubcats(self,orphanopt,writedir,zmin=-1.,zmax=1000.,ddir='/home/howdiedoo/boss/mksamplecatslatestdr12/',sampletag='cmass',runtag='dr12v4',cattypetag='Reid-targ',sysopt=1,fkpopt=0,NorSopt=2):
    """
    Copied from bootstrapclean/writesubcats200.py, pay attention to orphan option.
    following convention of catfitstotxt to output text files for xi computation.
    """

    ## copying from catfitstotxt; should have merged this instead of copying code.  Fix later!
    ## enforce sanity on catalog type and weighting scheme to protect myself from screw-ups!
    cpopt = 2 # default, for regular Reid catalogs.
    targcat = 0
    if re.search('targ',cattypetag):
      targcat = 1  # no redshifts!
      assert fkpopt == 0
      cpopt == 0

    elif re.search('ang',cattypetag):
      assert fkpopt == 0
      cpopt == 1

    if self.setup == 0:
      self.setuplookup()

    degtorad = np.pi/180.
    NSlist = ['N','S']

    assert NorSopt >= 0 and NorSopt <= 2
    if NorSopt == 2:
      NorSlist = [0,1]
    else:
      NorSlist = [NorSopt]  

    ## following convention in catfitstotxt
    catappend = ''
    if sysopt == 1:
      catappend = catappend + '-wsys'
    if fkpopt == 1:
      catappend = catappend + '-wfkp'

    ## keep account of some statistics of how well the boot-strap regions did.
    pixhists = np.zeros([2,self.nsub],dtype='int')
    tot = np.zeros(2,dtype='int')
    dropped = np.zeros(2,dtype='int')

    for DorR, fendin, fendout in zip([0,1],['.dat.fits','.ran.fits'],['.dat.txt','.ran.txt']):
      ofplist = []
      for ns in range(self.nsub):
        fout = writedir + sampletag + '-' + runtag + '-' + cattypetag + catappend + fendout+'.%04d' % (ns)
        ofplist.append(open(fout,'w'))

      for NorS in NorSlist:
        NStag = NSlist[NorS]
        fin = ddir + sampletag + '-' + runtag + '-' + NStag + '-' + cattypetag + fendin
        ## puts ra/dec into pixel coordinates (radians, shift by 90 degrees).
        ## z only filled in if it's in the original catalog, otherwise, filled with -1.
        ra,dec,z,wgt = readsortradecwgtfromfits(fin,DorR=DorR,decsortopt=0,cpopt=cpopt,sysopt=sysopt)
        ## make sure they're not sorted!  Need to output ra/dec in degrees.
        dd = fitsio.read(fin)
        assert (np.fabs(dec*180./np.pi - dd['DEC']) < 2.0e-6).all()

        for ii in range(len(ra)):
          mypix = self.getpix(NorS,ra[ii],dec[ii],orphanopt)
          if mypix == -1:
            dropped[DorR] += 1
          else:
            if mypix < 0 or mypix >= self.nsub:
              print 'wtf?!',mypix,NorS,ra[ii],dec[ii]
            assert mypix >= 0 and mypix < self.nsub
            pixhists[DorR][mypix] += 1
            if targcat == 0:
              if z[ii] < zmin or z[ii] > zmax:
                continue
              ofplist[mypix].write('%.12e %.12e %.6e %.6e\n' % (dd['RA'][ii],dd['DEC'][ii],dd['Z'][ii],wgt[ii]))
            else: #target catalog, no z!
              ofplist[mypix].write('%.12e %.12e %.6e\n' % (dd['RA'][ii],dd['DEC'][ii],wgt[ii]))

      for ns in range(self.nsub):
        (ofplist[ns]).close()

    print 'this many dropped: ',dropped
    print 'this many total: ',tot
    print 'frac of randoms kept:',float(dropped[1])/float(tot[1])
    ## stats on R hists:
    pixR = (pixhists[1,:]).copy()
    pixR.sort()
    rkept = pixR.sum()
    if not (rkept+dropped[1] == tot[1]):
      print 'beth, this did not work out.',rkept,dropped[1],tot[1]
    mymed = 0.5*(pixR[-self.nsub/2]+pixR[-self.nsub/2-1])
    print 'median occupancy',mymed
    print 'fractional var about median',(((pixR-mymed)**2).sum()/float(len(pixR)))**0.5/mymed



def plotsubcats(Nsub, fbase, ax=None):
  """
  Plot subregions by cycling over different colors to visualize them.
  fbase should be full path to the subregion files, minus '.%04d' at the end to specify each subregion.
  """
  if ax is None:
    ff = plt.figure(figsize=[6,6])
    ax=ff.add_subplot(1,1,1)
  else:
    ff = None

  clrlist = ['b.','g.','r.','c.','m.','y.']

  for ii in range(Nsub):
    #if ii > 10: break
    aa = np.loadtxt(fbase+'.%04d' % (ii),usecols=[0,1])
    ra = aa[:,0]
    ra = ra - 90.0
    xx = np.where(ra < 0.0)[0]
    ra[xx] += 360.0
    del xx
    clr = clrlist[ii%(len(clrlist))]
    ax.plot(ra,aa[:,1],clr)


class bootsec:  ## also divide survey into bootstrap regions generated by joining sectors together up to some area.
  def __init__(self):
    pass


