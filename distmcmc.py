import numpy as np
import os
import re
import scipy.interpolate as interp
import cosmo
import sys
import mcmcutils
import time
import scipy.optimize

### this is the reach of the sound horizon interpolator, don't integrate past that!
global och2min, och2max, doch2, noch2, obh2min, obh2max, dobh2, nobh2
global rsdat ## grid for zstar, rstar, rdrag interpolation.
global mCMB, icovCMB  ## mean and cov of 3x3 CMB distance priors.

global abaolist, DVorsfidlist ## for z =0.57, 0.32
global rsfid, DAorsfidlist, Hrsfidlist  ## for anisotropic.

global distparams # = ['omegabh2','omegach2','omegak','w','wa','H0']  ## list of cosmo parameters you will vary
global distdefaults
global daperp, dapar, aperpmin, aparmin, naperp, napar, anidat

def getcambtable(och2list,obh2list,onuh2list=[0.0006450616], nusplitopt = 0, ofp=None):
  """
  Runs camb at every value of och2, obh2 and grabs zstar, rs(zstar), zdrag, rs(zdrag).  
  Run this from the directory containing camb.
  Just adding option for zooming over neutrino masses.  We could do single mass or degenerate as test cases.
  nusplitopt = 0: single mass state.
  nusplitopt = 1: degenerate.
  these two are not yet written.
  nusplitopt = 2: normal hierarchy.
  nusplitopt = 3: inverted hierarchy.
  """

  nget = 4

  cfname = "cambtmp.ini"
  if len(onuh2list) == 1:
    cmbdat = np.zeros([len(obh2list), len(och2list), nget+2])
  else:
    cmbdat = np.zeros([len(obh2list), len(och2list), len(onuh2list), nget+3])
      

  ## names we want from camb run
  taglist = ['^zstar', '^r_s\(zstar\)', '^zdrag', '^r_s\(zdrag\)']

  for (ii, obh2) in zip(range(len(obh2list)), obh2list):
    for (jj, och2) in zip(range(len(och2list)), och2list):
      for (kk, onuh2) in zip(range(len(onuh2list)), onuh2list):
        cfp = open(cfname,'w')
        cfp.write('DEFAULT(bethdefaultmnu_getrsinterp_v2.ini)\n')
        cfp.write('ombh2 = %e\n' % (obh2))
        cfp.write('omch2 = %e\n' % (och2)) 
        cfp.write('omnuh2 = %e\n' % (onuh2)) 

        if nusplitopt == 0:
          pass
        elif nusplitopt == 1:
          m1 = cosmo.onuh2tomdegenerate(onuh2)
          m1frac = 1./3.
          m2frac = 1./3.
          m3frac = 1./3.
        elif nusplitopt == 2:
          m1tmp = cosmo.onuh2tomdegenerate(onuh2)
          if m1tmp*3 <= cosmo.Smnumin[0]:
            continue
          m1, m2, m3 = cosmo.setnumasses(m1tmp*3,0)
          mtot = m1 + m2 + m3
          m1frac = m1/mtot
          m2frac = m2/mtot
          m3frac = m3/mtot
          print 'hi beth Normal',m1,m2,m3,m1frac,m2frac,m3frac
        elif nusplitopt == 3:
          m1tmp = cosmo.onuh2tomdegenerate(onuh2)
          if m1tmp*3 <= cosmo.Smnumin[1]:
            continue
          m1, m2, m3 = cosmo.setnumasses(m1tmp*3,1)
          mtot = m1 + m2 + m3
          m1frac = m1/mtot
          m2frac = m2/mtot
          m3frac = m3/mtot
          print 'hi beth inverted',m1,m2,m3,m1frac,m2frac,m3frac
        else:
          print 'wacky nuopt!'
          sys.exit(1)
         
        if nusplitopt == 0:  ## one mass eigenstate.
          ## copied from default camb file.
          cfp.write('massless_neutrinos = 2.03066666667\n')
          cfp.write('massive_neutrinos = 1\n')
          cfp.write('nu_mass_eigenstates = 1\n')
          cfp.write('nu_mass_degeneracies = 1.01533333333\n')
          cfp.write('nu_mass_fractions = 1\n')
        else: # 3 mass eigenstates.
          cfp.write('massless_neutrinos = 0.0\n')
          cfp.write('massive_neutrinos = 3.046\n')
          cfp.write('nu_mass_eigenstates = 3\n')
          cfp.write('nu_mass_degeneracies = 1 1 1\n')
          cfp.write('nu_mass_fractions = %.6e %.6e %.6e\n' % (m1frac, m2frac, m3frac))

        cfp.close()
        print obh2, och2, onuh2
        # tmp run twice!!!
        #mystr = './camb %s' % (cfname)
        #print mystr
        #os.system(mystr)
        print 'take tmp first camb call out!'
        mystr = './camb %s > tmpout' % (cfname)
        os.system(mystr)

        ifp = open('tmpout','r')
        vallist = np.zeros(nget) - 1.
        for line in ifp:
          for ti in range(len(taglist)):
            if re.search(taglist[ti],line):
              assert vallist[ti] < -0.5
              vallist[ti] = float(line.strip('\n').split('=')[1])
        ifp.close()
        assert (vallist > 0.).all()
        if len(onuh2list) == 1:
          cmbdat[ii,jj,0] = obh2
          cmbdat[ii,jj,1] = och2
          cmbdat[ii,jj,2:] = vallist[:]
          if(ofp is not None):
            for qq in range(nget+2):
              ofp.write('%e ' % (cmbdat[ii,jj,qq]))
            ofp.write('\n')
        else:
          cmbdat[ii,jj,kk,0] = obh2
          cmbdat[ii,jj,kk,1] = och2
          cmbdat[ii,jj,kk,2] = onuh2
          cmbdat[ii,jj,kk,3:] = vallist[:]
          if(ofp is not None):
            for qq in range(nget+3):
              ofp.write('%e ' % (cmbdat[ii,jj,kk,qq]))
            ofp.write('\n')


def printcambtable(cmbdat,ofp,onuh2list=[0.0006450616]):
  """
  not actually used??
  """
  if len(onuh2list) == 1:
    (nx, ny, nz) = cmbdat.shape()
    for xi in nx:
      for yi in ny:
        for zi in nz:
          ofp.write('%e ' % (cmbdat[xi,yi,zi]))
        ofp.write('\n')
  else:
    (nw, nx, ny, nz) = cmbdat.shape()
    for wi in nw:
      for xi in nx:
        for yi in ny:
          for zi in nz:
            ofp.write('%e ' % (cmbdat[xi,yi,zi]))
          ofp.write('\n')


def anibaosetup(cmassanifname):
##Antonio, the file in the attachment contains the consensus constraints
#for DR11 post-reconstruction with the same format as the one you sent
#us, that is
#
#   alpha_perp, alpha_para, P(alpha_perp, alpha_para)
  
  global daperp, dapar, aperpmin, aparmin, naperp, napar, anidat

  try:
    dah = np.loadtxt(cmassanifname)
  except:
    print 'need to setup cmass ani in this directory, try again!'
    sys.exit(1)
  ncol = len(dah[0,:])

  naperp = len(np.where(dah[:,1] == dah[0,1])[0])
  napar = len(np.where(dah[:,0] == dah[0,0])[0])
  aperp1d = dah[:,0].reshape(naperp,napar)[:,0]
  apar1d = dah[:,1].reshape(naperp,napar)[0,:]
  daperp = (aperp1d[1:] - aperp1d[:-1]).mean()
  dapar = (apar1d[1:] - apar1d[:-1]).mean()
  aperpmin = aperp1d.min()
  aparmin = apar1d.min()

  assert (np.fabs(aperp1d[1:] - aperp1d[:-1] - daperp) < 2.0e-6).all()
  assert (np.fabs(apar1d[1:] - apar1d[:-1] - dapar) < 2.0e-6).all()
  print 'passed linear grid check for Ariel anisotropic.'

  ## sanity check on P -> chi2.  didn't forget a factor of 2 right??

  ## take out the 0.'s in the file, replace with a really big chi2.
  tt = np.where(dah[:,2] > 0.)[0]
  mymin = dah[tt,2].min()
  mymax = dah[tt,2].max()
  tt = np.where(dah[:,2] < 0.5*mymin)[0]
  qq = np.where(dah[:,2] >= 0.5*mymin)[0]
  dah[tt,2] = -2.*np.log(mymin) + 100. 
  dah[qq,2] = -2.*np.log(dah[qq,2])

  anidat = dah.reshape(naperp,napar,ncol)

  xx = np.where(aperp1d > 1.000)[0]
  yy = np.where(apar1d > 1.000)[0]

def rssetup():
  global och2min, och2max, doch2, noch2, obh2min, obh2max, dobh2, nobh2
  global rsdat
  try:
    rsdattmp = np.loadtxt('rscambfine.dat',usecols=[0,1,2,3,4,5])
  except:
    print 'need to set up rscambfine.dat in this directory, try again!'
    sys.exit(1)
  ncol = len(rsdattmp[0,:])
  noch2 = len(np.where(rsdattmp[:,0] == rsdattmp[0,0])[0])
  nobh2 = len(np.where(rsdattmp[:,1] == rsdattmp[0,1])[0])
    
  rsdat = rsdattmp.reshape(nobh2,noch2,ncol)


  ## set global bounds.
  ## check for linear spacing.

  och2min = rsdat[0,:,1].min()
  och2max = rsdat[0,:,1].max()
  doch2 = (rsdat[0,1:,1] - rsdat[0,:-1,1]).mean()
  noch2 = len(rsdat[0,:,1])
  chk = (rsdat[0,1:,1] - rsdat[0,:-1,1]).std()
  assert np.fabs(chk) < 2.0e-6
  obh2min = rsdat[:,0,0].min()
  obh2max = rsdat[:,0,0].max()
  dobh2 = (rsdat[1:,0,0] - rsdat[:-1,0,0]).mean()
  nobh2 = len(rsdat[1:,0,0])
  chk = (rsdat[1:,0,0] - rsdat[:-1,0,0]).std()
  assert np.fabs(chk) < 2.0e-6
  print obh2min, obh2max, dobh2, nobh2, och2min, och2max, doch2, nobh2

def bilinearani(aperptmp,apartmp):
  global daperp, dapar, aperpmin, aparmin, naperp, napar, anidat

  aperp = aperptmp
  apar = apartmp
  wantscalar = 0

  try:
    ll = len(aperp)
  except:
    aperp = np.array([aperp])
    apar = np.array([apar])
    wantscalar = 1
  if len(aperp) != len(apar):
    return -1.

  ix = np.array(np.floor((aperp - aperpmin)/daperp),dtype='int')
  iy = np.array(np.floor((apar - aparmin)/dapar),dtype='int')

  tmp = np.where(ix < 0)[0]
  ix[tmp] = 0
  tmp = np.where(ix > naperp-2)[0]
  ix[tmp] = naperp-2

  tmp = np.where(iy < 0)[0]
  iy[tmp] = 0
  tmp = np.where(iy > napar-2)[0]
  iy[tmp] = napar-2

  fx = (aperp - aperpmin)/daperp - ix
  fy = (apar - aparmin)/dapar - iy

  result = (1.-fx)*(1.-fy)*anidat[ix,iy,2] + \
           (1.-fx)*(fy)*anidat[ix,iy+1,2] + \
           (fx)*(1.-fy)*anidat[ix+1,iy,2] + \
           (fx)*(fy)*anidat[ix+1,iy+1,2]
  if wantscalar == 1:
    result = result[0]
  return result

## bilinear interpolators.
def bilinear(obh2tmp, och2tmp,whichI=0):
  """ whichI = 0 for rszstar, whichI = 1 for rzdrag, 2 for zstar, 3 for zdrag"""


  global och2min, och2max, doch2, noch2, obh2min, obh2max, dobh2, nobh2
  global rsdat

  obh2 = obh2tmp
  och2 = och2tmp
  wantscalar = 0
  try:
    ll = len(obh2)
  except:
    obh2 = np.array([obh2])
    och2 = np.array([och2])
    wantscalar = 1
  if len(obh2) != len(och2):
    return -1.

  ix = np.array(np.floor((obh2 - obh2min)/dobh2),dtype='int')
  iy = np.array(np.floor((och2 - och2min)/doch2),dtype='int')

  tmp = np.where(ix < 0)[0]
  ix[tmp] = 0
  tmp = np.where(ix > nobh2-2)[0]
  ix[tmp] = nobh2-2

  tmp = np.where(iy < 0)[0]
  iy[tmp] = 0
  tmp = np.where(iy > noch2-2)[0]
  iy[tmp] = noch2-2

  fx = (obh2 - obh2min)/dobh2 - ix
  fy = (och2 - och2min)/doch2 - iy

  if whichI == 0:
    iz = 3
  if whichI == 1:
    iz = 5
  if whichI == 2:
    iz = 2
  if whichI == 3:
    iz = 4
  
  result = (1.-fx)*(1.-fy)*rsdat[ix,iy,iz] + \
           (1.-fx)*(fy)*rsdat[ix,iy+1,iz] + \
           (fx)*(1.-fy)*rsdat[ix+1,iy,iz] + \
           (fx)*(fy)*rsdat[ix+1,iy+1,iz]
  if wantscalar == 1:
    result = result[0]
  return result

def rs(obh2tmp, och2tmp,whichrs=0):
  """ whichrs = 0 for rszstar, whichrs = 1 for rzdrag"""
  return bilinear(obh2tmp, och2tmp, whichI=whichrs)

def zCMB(obh2tmp, och2tmp, whichz=0):
  """ whichz = 0 for zstar, whichz = 1 for zdrag"""
  return bilinear(obh2tmp, och2tmp, whichI=whichz+2)

def dCMB(pcc):
  """
  Return vector of distance prior parameters R, la, obh2 from a cosmo instance.
  """
  myzstar = zCMB(pcc.obh2, pcc.och2,0)
  myDAzstar = pcc.DAz(1./(1.+myzstar))
  myrszstar = rs(pcc.obh2, pcc.och2,0)
  return np.array([np.sqrt(pcc.obh2 + pcc.och2)*myDAzstar/cosmo.DH, \
                   np.pi*myDAzstar/myrszstar,
                   pcc.obh2]), np.array([myzstar, myDAzstar, myrszstar])


def printPlanckprior(mvec,cov,outfname):
  ofp = open(outfname,'w')
  for i in range(3):
    ofp.write('%e\n' % (mvec[i]))
  for i in range(3):
    for j in range(3):
      ofp.write('%e ' % (cov[i,j]))
    ofp.write('\n')
  ofp.close()

def readPlanckprior(infname):
  ifp = open(infname,'r')
  mvec = np.zeros(3)
  i=0
  for line in ifp:
    mvec[i] = float(line.strip('\n'))
    i += 1
    if i == 3: break
  ifp.close()

  cov = np.loadtxt(infname,skiprows=3)
  return mvec, cov

def getPlanckprior(fbase,nchains=8,outfname=None,docheck=0):
  """
  3x3 Planck prior using parameters R, la, obh2 (see http://arxiv.org/pdf/1304.4514v2.pdf)
  But we're not using analytic expressions for zstar, zdrag, rs(z), etc.  We use CAMB outputs that
  depend only on och2 and obh2, and were computed here with minimal non-zero neutrino mass.
  """

  nchk = 500
  print 'Planck prior results for',fbase
  pfname = fbase+'.paramnames'
  mycnt = 0
  for ci in range(1,nchains+1):
    ftmp = fbase+'_%d.txt' % ci
    nn, cc = cosmo.readPlanckchain(pfname,ftmp) 
    myDAzstari = 100.*cc['rstar'][:]/cc['thetastar'][:]
    Ri = np.sqrt(cc['omegabh2'][:] + cc['omegach2'][:])*myDAzstari/cosmo.DH
    lai = np.pi*100./cc['thetastar'][:]
    #lai = cc['thetastar'] ## tmp!
    obh2i = cc['omegabh2'][:]
    wgti = cc['weight'][:]
    if ci == 1:
      R = Ri.copy()
      la = lai.copy()
      obh2 = obh2i.copy()
      wgt = wgti.copy()
      if docheck == 1: ## sanity check myDA vs theres.
        j=0
        mytmp = 0
        mydiff = np.zeros(3)
        while j < len(R):
          pcc = cosmo.Planck2cosmo(cc[j],nn)
          tCMB, extra = dCMB(pcc)
          mydiff += tCMB - np.array([R[j],la[j],obh2[j]])
          ### now returning extra information for plotting from the chains later.
          #mydiff += dCMB(pcc) - np.array([R[j],la[j],obh2[j]])
#          print j, dCMB(pcc), np.array([R[j],la[j],obh2[j]]), mydiff
          j += int(len(R)/float(nchk))
          mytmp += 1
        mydiff = mydiff/float(mytmp)
    else:
      R = np.concatenate((R,Ri))
      la = np.concatenate((la,lai))
      obh2 = np.concatenate((obh2,obh2i))
      wgt = np.concatenate((wgt,wgti))
    mycnt += len(Ri)

  assert len(R) == mycnt
  wgtsum = wgt.sum()
  mvec = np.zeros(3)
  cov = np.zeros([3,3])
  for pi, pp in zip(range(3), [R, la, obh2]):
    mvec[pi] = (pp*wgt).sum()/wgtsum
  print 'mean',mvec
  for pi, pp in zip(range(3), [R, la, obh2]):
    for qi, qq in zip(range(3), [R, la, obh2]):
      cov[pi,qi] = ((pp-mvec[pi])*(qq-mvec[qi])*wgt).sum()/wgtsum
  print 'err',[(cov[pi,pi])**0.5 for pi in range(3)]
  print 'normcov'
  for pi, pp in zip(range(3), [R, la, obh2]):
    for qi, qq in zip(range(3), [R, la, obh2]):        
      print cov[pi,qi]/(cov[pi,pi]*cov[qi,qi])**0.5

  if(docheck == 1):
    print 'mean fractional diffs:',mydiff/mvec
  else:
    print 'passed check already'
    print 'fractional accuracy on CMB calc: [  4.49937390e-05   4.41190131e-05   0.00000000e+00]'
  if outfname is not None:
    printPlanckprior(mvec,cov,outfname)
  return mvec, cov

def getPlanckstepmat(fparams, fbase, paramlist, outfname, steprescale=2.4,nchains=8):
  """
  return a step matrix for the input chain using mcmcutils.
  """
  ## copying from readPlanckchain, need a hack to make readable to mcmcutils.
  names = ['weight','lnlike']
  ifpp = open(fparams,'r')
  for line in ifpp:
    nn = line.split('\t')[0].strip(' ').strip('*')
    names.append(nn)
  ## write to a file
  ofp = open('planckcolstmp.dat','w')
  for i, nn in zip(range(len(names)), names):
    if i==0:
      ofp.write('%s' % (nn))
    else:
      ofp.write(',%s' % (nn))
  ofp.close()
  mystep = []
  for i in range(1,nchains+1):
    fchain = fbase + '_%d.txt' % i
    cc = mcmcutils.chain(fchain,'planckcolstmp.dat')
    ## hack -- don't vary anything except the parameters we care about.
    cc.mcmcfixed[:] = 1 ## set everything to fixed.
    for pp in paramlist:
      cc.mcmcfixed[cc.mcmcp[pp]] = 0

    cc.fillstepmatrix(steprescale=steprescale)
    if i==1:
      mystep = cc.step_mat.copy()
    else:
      mystep = mystep + cc.step_mat
  mystep = mystep/float(nchains)
  ## hack!!
  ## take the average and then print.
  cc.step_mat = mystep
#  cc.printstepmatrix(outfname)
  ofp = open(outfname,'w')
  ## print a list o the parameters.
  xx = np.where(cc.mcmcfixed == 0)[0]
#  print cc.mcmcpreverse[xx]
  orderednames = [cc.mcmcpreverse[xx[ii]] for ii in range(len(xx))]
  assert len(orderednames) == len(paramlist)
  for i in range(len(paramlist)):
    if i == 0:
      ofp.write('# %s' % (orderednames[i]))
    else:
      ofp.write(',%s' % (orderednames[i]))
  ofp.write('\n')
  for i in range(len(paramlist)):
    for j in range(len(paramlist)):
      ofp.write('%e ' % (mystep[i,j]))
    ofp.write('\n')

  ofp.close()

def readPlanckstepmat(infname):
  ifp = open(infname,'r')
  line = ifp.readline()
  ifp.close()
  nn = [ii.strip('# \n') for ii in line.split(',')]
  stepmat = np.loadtxt(infname,skiprows=1)
  return stepmat, nn


def CMBchi2(cc):
  """
  Input cosmology.  Output likelihood using CMB distance prior information.
  """
  global mCMB, icovCMB
  tCMB, extra = dCMB(cc)
  dvec = np.matrix(tCMB - mCMB)
  chi2 = ((dvec*icovCMB) * (dvec.T))[0,0]
  return chi2, extra

def BAOchi2(cc,whichbao):
  """
  Returns total BAO likelihood given a mask of which BAO likes to compute.
  whichBAO = [CMASS iso, LOWZ iso, CMASS ani] for now.
  """
  global DVorsfidlist ## for z =0.57, 0.32
  global abaolist
  ## 11/1 just added this line.
  global DAorsfidlist, Hrsfidlist


#  if whichbao[0] == 0 and whichbao[1] == 0 and whichbao[2] == 0:
### hope, this setup is overridden with hacked version of CMASSiso
#    return 0.

  ## hard code for now!!
  DVorslist = np.array([1.0144*DVorsfidlist[0], 1.018*DVorsfidlist[1]])
  ### LOWZ error went up by 0.001..
  #DVorssig2list = np.array([(0.0098*DVorsfidlist[0])**2, (0.020*DVorsfidlist[1])**2])
  DVorssig2list = np.array([(0.0098*DVorsfidlist[0])**2, (0.021*DVorsfidlist[1])**2])


  ## hack wrong alpha center.  
  if whichbao[0] == 0 and whichbao[1] == 0 and whichbao[2] == 0:
    DVorslist = np.array([1.027*DVorsfidlist[0], 1.018*DVorsfidlist[1]])

  ## compute all desired quantities up front, put into baoinfo.
  rsval = rs(cc.obh2, cc.och2, 1)
  BAOinfo = np.array([[cc.DVMpc(abaolist[0])/rsval, cc.DAz(abaolist[0])/rsval, cc.Hofz(abaolist[0])*rsval], \
                      [cc.DVMpc(abaolist[1])/rsval, cc.DAz(abaolist[1])/rsval, cc.Hofz(abaolist[1])*rsval]])
         
  baochi2 = 0.
  if whichbao[0] == 0 and whichbao[1] == 0 and whichbao[2] == 0:
    i = 0
    baochi2 += ((BAOinfo[i,0] - DVorslist[i])**2/DVorssig2list[i])

  ## isotropic.
  for i in range(2):
    if whichbao[i] == 1:
      baochi2 += ((BAOinfo[i,0] - DVorslist[i])**2/DVorssig2list[i])

  ## anisotropic:
  if whichbao[2] == 1:
    ## convert model DA and H to aperp and apar.
    aperp = BAOinfo[0,1]/DAorsfidlist[0]
    apar = Hrsfidlist[0]/BAOinfo[0,2]
    chi2ani = bilinearani(aperp,apar)
    baochi2 += chi2ani

  ## hack!!
  if whichbao[0] == 3 and whichbao[1] == 3 and whichbao[2] == 3:
    ## map (1.027, 1.027) to the minimum chi2 point.
    aperpmin = 1.04571428
    aparmin = 0.96485717
    aperp = BAOinfo[0,1]/DAorsfidlist[0] - 1.027 + aperpmin
    apar = Hrsfidlist[0]/BAOinfo[0,2] - 1.027 + aparmin
    chi2ani = bilinearani(aperp,apar)
    baochi2 += chi2ani
    


  return baochi2, np.concatenate(BAOinfo.flatten(), np.array([rsval]))

def setuplikelihoods(cmbfname,cmassanifname=None):
  global mCMB, icovCMB
  global abaolist, DVorsfidlist
  global DAorsfidlist, Hrsfidlist
  global rsfid

  abaolist = np.array([1./1.57, 1./1.32])
  mCMB, covtmp= readPlanckprior(cmbfname)
  icovCMB = (np.matrix(covtmp)).I
  DVorsfidlist = np.array([2])
  obh2fid = 0.0224
  och2fid = 0.11186
  hfid = 0.7
  ccfid = cosmo.cosmo(och2=och2fid,obh2=obh2fid,h=hfid,forceflat=1)
  rsfid = rs(obh2fid,och2fid,1)
  DV1fid = ccfid.DVMpc(abaolist[0])
  DV2fid = ccfid.DVMpc(abaolist[1])

  ## set up anisotropic fiducial values.
  DA1fid = ccfid.DAz(abaolist[0])
  DA2fid = ccfid.DAz(abaolist[1])
  H1fid = ccfid.Hofz(abaolist[0])
  H2fid = ccfid.Hofz(abaolist[1])

  DVorsfidlist = np.array([DV1fid/rsfid, DV2fid/rsfid])
  DAorsfidlist = np.array([DA1fid/rsfid, DA2fid/rsfid])
  Hrsfidlist = np.array([H1fid*rsfid, H2fid*rsfid])


  print 'BAO stuff:',DV1fid, DV2fid, rsfid
  print DVorsfidlist
  print DAorsfidlist
  print Hrsfidlist
  print 'chi2 for fiducial cosmology'
  cmbchi2val, extra = CMBchi2(ccfid)
  baochi2val, extrabao = BAOchi2(ccfid,np.array([1,1,0],dtype='int'))

  print 'CMB: ',cmbchi2val
  print 'BAO: ',baochi2val

  if cmassanifname is not None:
    anibaosetup(cmassanifname)
    baochi2valani, extrabao = BAOchi2(ccfid,np.array([0,0,1],dtype='int'))
  print 'BAOani: ',baochi2valani


## copying from bethalexie/mcmc.c

def mcmcstep(old, stepmat):

  #tmp = (np.random.normal(size=10000))
  ### works!
  #print 'chk gauss',tmp.mean(), tmp.std()
  step = np.matrix(np.random.normal(size=len(old)))
  sig = stepmat*step.T
  new = np.array(old + sig.T)
  new = np.array(new[0,:])
  return new
  
def chain2cosmo(elt, eltnames, mcmcfixed):
## copying Planck2cosmo in cosmo.py
  global distdefaults
  global distparams

## using defaults in distdefaults
  mydist = distdefaults.copy()
  j = 0
  for i in range(len(distparams)):
    if mcmcfixed[i] == 0:  ## then set the value from elt.
      mydist[distparams[i]] = elt[j]
      j += 1
  assert j == len(elt)

  h = mydist['H0']*0.01
  och2 = mydist['omegach2']
  obh2 = mydist['omegabh2']
  okh2 = mydist['omegak']*h**2
  w0 = mydist['w']
  wa = mydist['wa']

  ogh2 = 2.469e-5  ## default value givne Tcmb.
  onuh2val = 0.0006450616  ## taken from camb ini file Mar13 base_planck_lowl_lowLike.ini; this is z=0 value.
  SorD = 0 # single mass eigenstate.

  odeh2 = h**2 - (och2 + obh2 + ogh2 + onuh2val + okh2)
  oDE = odeh2/h**2

  if 'omegak' in eltnames:
    forceflat = 0
  else:
    forceflat = 1

  ## set neutrino defaults:
  cc = cosmo.cosmo(och2=och2,\
             obh2=obh2,
             h=h,\
             w0=w0,\
             wa=wa,\
             omegaDE=oDE,\
             forceflat=forceflat,\
             onuh2val=onuh2val,\
             SorD=SorD)

  ## check omegak
  assert np.fabs(cc.ok - mydist['omegak']) < 2.0e-6
  return cc

def stdsetup():
  rssetup()
#  setuplikelihoods(cmbfname='base_planck_lowl_lowLike_highL.3x3',cmassanifname='dah_consensus_dr11_rec.dat')
  setuplikelihoods(cmbfname='base_planck_lowl_lowLike_highL.3x3',cmassanifname='dah_consensus_dr11_rec_sysfinal.dat')



def runchain(stepfname,mcmcfixed,nmax,chainfname,restartopt=0,whichbao=[1,0],distdefaultsinput=None):
  """
  mcmcfixed should be length of distparams and specify which are varying.
  This should agree with step matrix -- we'll check that inside.
  If distdefaults = None, they will be set here according to my will.
  """
  global distparams
  global distdefaults

  ## set up stuff.
  stdsetup()

  ## fixed global list of cosmo parameters you will vary
  distparams = ['omegabh2','omegach2','omegak','w','wa','H0'] 
  ## for physical matter densities, take from
  ## base_planck_lowl_lowLike_highL_planckgausslikelihood.dat
  if distdefaultsinput is None:
    distdefaults = {'w':-1., 'wa':0.,'omegabh2':2.20745e-02,'omegach2':1.19802e-01,'omegak':0, 'H0':70.0}
  else:
    distdefaults = distdefaultsinput.copy()

  stepmat, eltnames = readPlanckstepmat(stepfname)
  nparam = len(stepmat[:,0])
  distindx = np.zeros(nparam,dtype='int') ## make sure distparams and stepmat have parameters in the same order.
  stepmat = np.matrix(stepmat)
  ## check
  si = 0
  for pi, p in zip(range(len(distparams)), distparams):
    if p in eltnames:
      assert mcmcfixed[pi] == 0
      distindx[si] = pi
      si += 1
    else:
      assert mcmcfixed[pi] == 1
  assert si == len(distindx)
  assert (distindx[:-1] < distindx[1:]).all()
  print distindx
  print 'passed mcmcfixed checks' 

  old = np.zeros(nparam)
  new = np.zeros(nparam)

  if restartopt == 1:
    print 'not written yet'
    sys.exit(1)

  else:
    ## fill in matrix that's varying.
    ## initialize with defaults
    oldd = distdefaults.copy() ## same size as distparams.
    for i in range(nparam):
      old[i] = oldd[eltnames[i]]

    oldcc = chain2cosmo(old,eltnames,mcmcfixed)
#def chain2cosmo(elt, eltnames, mcmcfixed):
    ## get old chi2.
    
    oldbaochi2, oldbaoinfo = BAOchi2(oldcc, whichbao)
    oldcmbchi2, oldcmbinfo = CMBchi2(oldcc)
    oldchi2 = oldbaochi2 + oldcmbchi2

    ## print some info at the top of the chain
    cfp = open(chainfname,'w')
    cfp.write('# stepfile: %s\n' % stepfname)
    for i in range(len(eltnames)):
      if i == 0:
        cfp.write('# %s' % (eltnames[i]))
      else:
        cfp.write(',%s' % (eltnames[i]))
    cfp.write('\n')
    cfp.close()

  print 'starting cosmo',old,oldbaochi2,oldcmbchi2,oldchi2
  print oldcc

  naccept = 0  
  nreject = 0
  tstart = time.time()
  currwgt = 0
  while naccept < nmax:  
    new = mcmcstep(old,stepmat)
    newcc = chain2cosmo(new,eltnames,mcmcfixed)
    ## get new chi2.
    newbaochi2, newbaoinfo = BAOchi2(newcc, whichbao)
    newcmbchi2, newcmbinfo = CMBchi2(newcc)
    newchi2 = newbaochi2 + newcmbchi2

#    print 'newcosmo',new,newbaochi2,newcmbchi2,newchi2
#    print newcc
    rr = np.random.random()
    qq = np.exp(-0.5*(newchi2 - oldchi2))
    #if(newchi2 <= oldchi2 or (np.random.random() < np.exp(-0.5*(newchi2 - oldchi2)))):
    if(newchi2 <= oldchi2 or (rr < qq)):
      ## open and close file pointer so we don't lose anything along the way, can inspect file
      ## while running!

#      printChainElement(old,cfp)
      cfp = open(chainfname,'a')
      cfp.write('%d %e %e %e ' % (currwgt,oldbaochi2,oldcmbchi2,oldchi2))
      for i in range(nparam):
        cfp.write('%e ' % (old[i]))
      for i in range(len(oldcmbinfo)):
        cfp.write('%e ' % (oldcmbinfo[i]))
      for i in range(len(oldbaoinfo.flatten())):
        cfp.write('%e ' % (oldbaoinfo.flatten()[i]))
      cfp.write('\n')
      cfp.close()
      old = new
      oldcc = newcc
      oldbaochi2 = newbaochi2
      oldcmbchi2 = newcmbchi2
      oldchi2 = newchi2
      oldcmbinfo = newcmbinfo.copy()
      oldbaoinfo = newbaoinfo.copy()

      naccept += 1
      currwgt = 0
    else:
      currwgt += 1
      nreject += 1
    if naccept % 100 == 0:
      t2 = time.time()
      print 'accepted %d, rejected %d, this took %e seconds.\n' % (naccept,nreject,t2-tstart) 

## print last chain elt.
#printChainElement(old,cfp)
  cfp = open(chainfname,'a')
  cfp.write('%d %e %e %e ' % (currwgt,oldbaochi2,oldcmbchi2,oldchi2))
  for i in range(nparam):
    cfp.write('%e ' % (old[i]))
  for i in range(len(oldcmbinfo)):
    cfp.write('%e ' % (oldcmbinfo[i]))
  for i in range(len(oldbaoinfo.flatten())):
    cfp.write('%e ' % (oldbaoinfo.flatten()[i]))
  cfp.write('\n')
  cfp.close()


def mysolvela(h,cozold, \
    Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
    w0=-1.,wa=0.,omegak=0.):
  """
  copying from mysolveDAcurved.
  """
  if np.fabs(omegak) < 2.0e-6:
    forceflat = 1
  else:
    forceflat = 0
  okh2 = omegak*h**2
  oDEh2 = h**2 - (cozold.och2 + cozold.obh2 + cozold.ogh2 + cozold.onuh2(1.) + okh2)
  oDE = oDEh2/h**2
  coztmp = cosmo.cosmo(och2=cozold.och2, obh2=cozold.obh2, ogh2=cozold.ogh2, Tcmb=cozold.Tcmb, \
                 h=h,
                 w0=w0,
                 wa=wa,
                 omegaDE = oDE,
                 forceflat=forceflat,
                 Neff=Neff,mnu1=mnu1,mnu2=mnu2,mnu3=mnu3,mnu4=mnu4,NorI=NorI,
                 Smnu=Smnu,onuh2val=onuh2val,SorD=SorD)
  tCMB, extra = dCMB(coztmp)
  return tCMB[1] ## this is la


def cosmofromlacurved(cozold,latarget,\
        Neff=3.046,mnu1=0.,mnu2=0.,mnu3=0.,mnu4=0.,NorI=-1,Smnu=-1,onuh2val=-1.,SorD=1,\
        w0=-1.,wa=0.,omegak=0.):
  try:
    newhval = scipy.optimize.brentq(lambda x: mysolvela(x,cozold, Neff, mnu1, mnu2, mnu3, mnu4, NorI, Smnu, onuh2val,SorD,w0,wa,omegak) - latarget, cozold.h*0.5,cozold.h*2.)
  except:
    print 'new h value not found in factor of 2 of old h value.  recode or error!'
    return None

  ## success!
  if np.fabs(omegak) < 2.0e-6:
    forceflat = 1
  else:
    forceflat = 0
  okh2 = omegak*newhval**2
  oDEh2 = newhval**2 - (cozold.och2 + cozold.obh2 + cozold.ogh2 + cozold.onuh2(1.) + okh2)
  oDE = oDEh2/newhval**2

  coznew = cosmo.cosmo(och2=cozold.och2, obh2=cozold.obh2, ogh2=cozold.ogh2, Tcmb=cozold.Tcmb, \
                 h=newhval,
                 w0=w0,
                 wa=wa,
                 omegaDE = oDE,
                 forceflat=forceflat,
                 Neff=Neff,mnu1=mnu1,mnu2=mnu2,mnu3=mnu3,mnu4=mnu4,NorI=NorI,
                 Smnu=Smnu,onuh2val=onuh2val,SorD=SorD)
  tCMB, extra = dCMB(coznew)
  assert((tCMB[1] - latarget)/latarget < 1.0e-5)
  return coznew

## Daniel exercise --> what range of och2 are implied at fixed l_A, och2?
def pBAOoch2(cozold,zlist):
  """
  Outputs are not normalized, so only relative probabilities are correct.
  Input cosmology sets the defaults.
  returns all the alphas at list of input zvals
  NOT WRITTEN YET.
  """
  coznew = cosmofromlacurved(cozold,mCMB[1])
  return 0 



if __name__ == '__main__':


  if(len(sys.argv) != 4 and len(sys.argv) != 2):
    print 'Usage: python distmcmc.py model baoset nmcmc'
    print 'Usage: python distmcmc.py nusplitopt [not used, hard-coded to do all for now.]'
    print 'model [0-4]:'
    print 'LCDM,oLCDM,wCDM,w0waCDM,ow0wacdm,wCDM w fixed ocbh2'
    print 'baoset [0-4]:'
    print '[CMASS iso, LOWZ iso, CMASS ani]'
    print '0: [1,0,0]' 
    print '1: [0,1,0]' 
    print '2: [1,1,0]' 
    print '3: [0,0,1]' 
    print '4: [0,1,1]' 
    print 'nmcmc = number of chain elements (not trials) you want'
    sys.exit(1)

  if len(sys.argv) == 2:

    ff = '/home/howdiedoo/Planck/PLA/base/planck_lowl/base_planck_lowl_planckgausslikelihood.dat'
    ifp = open(ff,'r')
    line = ifp.readline()  ##nparam
    line = ifp.readline()
    och2cen = float(line.strip())
    line = ifp.readline()
    obh2cen = float(line.strip())
    line = ifp.readline()
    nscen = float(line.strip())
    icov = np.loadtxt(ff,skiprows=4)
    icov = np.matrix(icov)
    cov = icov.I
    och2sig = (cov[0,0])**0.5
    obh2sig = (cov[1,1])**0.5
    r = cov[0,1]/(och2sig*obh2sig)
    print och2sig, obh2sig, r
    print och2sig/och2cen
    print obh2sig/obh2cen

    ## run fine grid later.
    doch2 = 0.1*och2sig
    dobh2 = 0.1*obh2sig

    nsig = 2.

    ## run coarse grid to get started.
    doch2 = 0.5*och2sig
    dobh2 = 0.5*obh2sig

    och2list = np.arange(och2cen - nsig*och2sig, och2cen + (nsig + 0.5*doch2)*och2sig, doch2)
    obh2list = np.arange(obh2cen - nsig*obh2sig, obh2cen + (nsig + 0.5*dobh2)*obh2sig, dobh2)

    ## tmp!!
#    och2list = np.array([och2cen])
#    obh2list = np.array([obh2cen])

    onuh2list = np.arange(0.0006451,0.0111,0.0006451)
    print onuh2list
#    onuh2list = np.array([0.002,0.003])

    for tag, nusplitopt in zip(['single','degenerate','normal','inverted'], [0,1,2,3]):
      ofp = open('rscambmutmp_%s.dat' % tag,'w')
      ofp.write('# obh2, och2, onuh2, zstar, rs(zstar), zdrag, rs(zdrag)\n')
      getcambtable(och2list,obh2list,onuh2list=onuh2list, nusplitopt = nusplitopt, ofp=ofp)
      ofp.close()
    print 'finished running camb with neutrinos, exiting'
    sys.exit(1)


  ## this is for running quickMCMC.
  whichmodel = int(sys.argv[1])
  assert whichmodel >= 0 and whichmodel <= 5
  whichbaoset = int(sys.argv[2])
  nmax = int(sys.argv[3])
  assert whichbaoset >= 0 and whichbaoset <= 6
  if whichbaoset == 0: whichbao = np.array([1,0,0]); baotag = 'CMASSiso'
  if whichbaoset == 1: whichbao = np.array([0,1,0]); baotag = 'LOWZiso'
  if whichbaoset == 2: whichbao = np.array([1,1,0]); baotag = 'CMASSiso_LOWZiso'
  if whichbaoset == 3: whichbao = np.array([0,0,1]); baotag = 'CMASSani'
  if whichbaoset == 4: whichbao = np.array([0,1,1]); baotag = 'LOWZiso_CMASSani'
  if whichbaoset == 5: whichbao = np.array([0,0,0]); baotag = 'CMASSisohack'
  if whichbaoset == 6: whichbao = np.array([3,3,3]); baotag = 'CMASSanihack'

  #distparams = ['omegabh2','omegach2','omegak','w','wa','H0'] 
  if whichmodel == 0:
    mbase = 'base'
    mcmcfixed=np.array([0,0,1,1,1,0])
  if whichmodel == 1:
    mbase = 'base_omegak'
    mcmcfixed=np.array([0,0,0,1,1,0])
  if whichmodel == 2:
    mbase = 'base_w'
    mcmcfixed=np.array([0,0,1,0,1,0])
  if whichmodel == 3:
    mbase = 'base_w_wa'
    mcmcfixed=np.array([0,0,1,0,0,0])
  if whichmodel == 4:
    mbase = 'base_omegak_w_wa'
    mcmcfixed=np.array([0,0,0,0,0,0])
    print 'no step for this one yet'
    sys.exit(1)
  if whichmodel == 5:
    mbase = 'base_w_fixobch2'
    mcmcfixed=np.array([1,1,1,0,1,0])

  stepfname='%s.step' % mbase
  chainfname = 'distmcmcchains/%s_%s.chain' % (mbase,baotag)  

  ###### this is code run in cambMar13/camb directory to get rs interpolators.
  ## from 
  ## /home/howdiedoo/Planck/PLA/base/planck_lowl/base_planck_lowl_planckgausslikelihood.dat
  if 0==1:
    ff = '/home/howdiedoo/Planck/PLA/base/planck_lowl/base_planck_lowl_planckgausslikelihood.dat'
    ifp = open(ff,'r')
    line = ifp.readline()  ##nparam
    line = ifp.readline()
    och2cen = float(line.strip())
    line = ifp.readline()
    obh2cen = float(line.strip())
    line = ifp.readline()
    nscen = float(line.strip())
    icov = np.loadtxt(ff,skiprows=4)
    icov = np.matrix(icov)
    cov = icov.I
    och2sig = (cov[0,0])**0.5
    obh2sig = (cov[1,1])**0.5
    r = cov[0,1]/(och2sig*obh2sig)
    print och2sig, obh2sig, r
    print och2sig/och2cen
    print obh2sig/obh2cen
  
    ## run fine grid later.
    doch2 = 0.1*och2sig
    dobh2 = 0.1*obh2sig
  
    nsig = 4.
  
    ## run coarse grid to get started.
    doch2 = 0.5*och2sig
    dobh2 = 0.5*obh2sig
    
  
    och2list = np.arange(och2cen - nsig*och2sig, och2cen + (nsig + 0.5*doch2)*och2sig, doch2)
    obh2list = np.arange(obh2cen - nsig*obh2sig, obh2cen + (nsig + 0.5*dobh2)*obh2sig, dobh2)
    print len(och2list), len(obh2list)
    ofp = open('rscamb.dat','w')
    ofp.write('# obh2, och2, zstar, rs(zstar), zdrag, rs(zdrag)\n') 
    getcambtable(och2list, obh2list, ofp)
    ofp.close()

  if 0==1:
    fbase = '/Users/bareid/work/montserratdata/Planck/PLA/base/planck_lowl_lowLike_highL/base_planck_lowl_lowLike_highL'
    outfname = fbase.split('/')[-1] + '.3x3'
    rssetup()
    getPlanckprior(fbase,nchains=8,outfname=outfname,docheck=0)
    #setuplikelihoods(cmbfname=outfname,cmassanifname='dah_consensus_dr11_rec.dat')
    setuplikelihoods(cmbfname=outfname,cmassanifname='dah_consensus_dr11_rec_sysfinal.dat')

  if 1==1:
    print 'do standard startup'
    stdsetup()
    ### stdsetup is equivalent to this.
    #rssetup()
    #setuplikelihoods(cmbfname='base_planck_lowl_lowLike_highL.3x3',cmassanifname='dah_consensus_dr11_rec.dat')

    #distparams = ['omegabh2','omegach2','omegak','w','wa','H0'] 
    #runchain(stepfname='base.step',mcmcfixed=np.array([0,0,1,1,1,0]),nmax=50000,chainfname="distmcmcchains/base_bao1_0.chain",restartopt=0,whichbao=[1,0],distdefaultsinput=None)
    runchain(stepfname=stepfname,mcmcfixed=mcmcfixed,nmax=nmax,chainfname=chainfname,restartopt=0,whichbao=whichbao)


