import numpy as np
import matplotlib.pyplot as plt
import sys
#import copy
#import scipy.optimize
import os
#import re
#from scipy.integrate import quad
#import scipy.interpolate as interp
import re
import mcmcutils


def runmodel(Mmin=1.1460142e+13,M1=1.57663077e+14,alpha=1.298719e+00,\
             Mcut=3.238620e+11,sigmalogM=3.40466e-01,\
             hvscale = 1.0, ihvscale = 1.0, cenvfrac = 0.0,\
             whichbox=1,cenopt=2,satopt=2, force_cenforsat=1,\
             ftag='tmprun',pbase='./',maskedopt = 0,
             whichfit = [0,1,0],runslow=0,COMVopt=0,
             aperp = 1.0,apar = 1.0, ## new!!
             use5002or5003=0,useprecp2=None,whichPB=0,writecat=None):  ## make 5002 the default, so old code still runs. 
  """
  pbase is the path to the directory where you want to run ./runall 
  whichbox = 0 for L0, 1 for MWhires
  whichfit is 0 (no) 1 (yes) for [lensing, xiell, wp]
  use5002or5003=0 for old version, 1 for new version (5003)
  for 5003, path should be to bethalexieprecpinterpNOFOFxiellwpwcov probably.
  useprecp2 = None if no precp2 available, useprecp2 = FASTP_BASE of precp2 case.
  whichPB = 0 is default, only one we have precompute counts for right now.
  writecat = catfname allows you to write out a catalog to catfname; otherwise none written.
  March 7: added aperp/apar as parameters.
  """

  ## change directories
  mycwd = os.getcwd() ## go back here at the end.
  os.chdir(pbase)

  parambase = 'pauto_base_newv2.bat'
  hodfname = 'hodparams/' + ftag + '.hod'
  hodnopath = ftag + '.hod'
  batfname = ftag + '.bat'
  outfiletag = ftag + 'HV%.3f_IHV%.3f_CENV%.3f' % (hvscale,ihvscale,cenvfrac)

  ofp = open(hodfname,'w')
  ofp.write("# comments here.\n")
  for pp in [Mmin, sigmalogM, M1, alpha, Mcut]:
    ofp.write('%e\n' % (pp))
  ofp.close()

  ifp = open(parambase,'r')
  ofp = open(batfname,'w')

  if runslow == 0:
    writecat = None ## can't write a catalog without slow.

  ofp.write('ALPHAPERP       %.6f\n' % (aperp))
  ofp.write('ALPHAPAR       %.6f\n' % (apar))

  if writecat is None:
    ofp.write('WRITE_CAT      0\n')
    ofp.write('CATFNAME       blah.cat\n')
  else:
    ofp.write('WRITE_CAT      1\n')
    ofp.write('CATFNAME       %s\n' % (writecat))

  if whichbox == 2:
    assert whichfit[0] == 0 ## lensing
    assert use5002or5003 == 1
    assert whichPB >= 0 and whichPB <= 2
    if COMVopt > 0 and useprecp2 is None: print 'COMVopt not supported yet. Break!'; return None
    if maskedopt == 0 and runslow != 1: print 'maskedopt = 1 only supported.  Break!'; return None

  assert whichfit[0] == 0 ## lensing
  assert whichbox != 0  ## don't think this works, need to set up halo catalog and mass fxn properly through precompute.

  if whichbox == 0 or whichbox == 1:
    assert whichbox != 0
    ofp.write('FASTP_FBASE_WP   precomputeMWhiresmaskedcombo_HV1.000_IHV1.000_CENV0.000\n') # not used.
  if whichbox == 2:
    ofp.write('FASTP_FBASE_WP        precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n' % (whichPB))

  if (use5002or5003==0 and whichfit[2] == 1) or (use5002or5003==1 and whichfit[2] == 1 and whichfit[1] == 0):  ## wp
    #if whichbox == 0 or whichbox == 1:  ## whichbox = 0 doesn't use this.  wrong! uses a mass function.
    if whichbox == 1:  
      ofp.write('FASTP_FBASE           precompute_MWhiresplanck\n')
## moved up.
#      ofp.write('FASTP_FBASE_WP   precomputeMWhiresmaskedcombo_HV1.000_IHV1.000_CENV0.000\n') # not used.
    if whichbox == 2:
      ofp.write('FASTP_FBASE           precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n' (whichPB))
#      ofp.write('FASTP_FBASE_WP        precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n' % (whichPB))
    if runslow == 1:
      ofp.write('USE_FASTP_COUNTS      0\n')
    else:
      ofp.write('USE_FASTP_COUNTS      2\n')
    ## these wont be used anyways.
    if(use5002or5003 == 0):
      ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5002.txt\n')
      ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1_splits5_6.syswtheory\n')
      ofp.write('XIELL_REBINFNAME    bin1.txt.reformat.new\n')
      ofp.write('ELLMAX              2\n')
    else:
      ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5003bin2.txt.dummy\n')
      ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1fineMU_splits1_1.syswtheory.dummy\n')
      ofp.write('XIELL_REBINFNAME    bin1fineMU.txt.reformat.new\n')
      ofp.write('ELLMAX              2\n')


  else: ## do xiell
    if whichbox == 1:
      if maskedopt == 0:
        if COMVopt == 0:
          ## in fixed HV,IHV case, just use 2d grid.
          if (np.fabs(hvscale - 1.0) < 1e-4 and \
             np.fabs(ihvscale - 1.0) < 1e-4 and \
             np.fabs(cenvfrac) < 1e-4) or (runslow == 1):
             ofp.write('FASTP_FBASE    precomputeMWhirescombo_HV1.000_IHV1.000_CENV0.000\n')
             if runslow == 1:
               ofp.write('USE_FASTP_COUNTS      0\n')
             else:
               ofp.write('USE_FASTP_COUNTS    2\n')
          elif useprecp2 is not None:
            ofp.write('FASTP_FBASE        %s\n' % (useprecp2))
            ofp.write('USE_FASTP_COUNTS    2\n')
  
          else:  ## use the grid.
            ofp.write('FASTP_FBASE         gridv0/precomputeMWhirescombo.precptable\n')
            assert runslow == 0
  #          if runslow == 1:
  #            ofp.write('USE_FASTP_COUNTS      0\n')
  #          else:
            ofp.write('USE_FASTP_COUNTS    3\n')
  
        else:  ##COMVopt == 1
          #### add this option later?
          #elif useprecp2 is not None:
          #  ofp.write('FASTP_FBASE        %s\n' % (useprecp2))
          #  ofp.write('USE_FASTP_COUNTS    2\n')
          assert np.fabs(cenvfrac) < 0.001 or np.fabs(cenvfrac - 0.3) < 0.001 or runslow == 1
          if np.fabs(cenvfrac - 0.3) < 0.001:
            ofp.write('FASTP_FBASE       precomputeMWhiresCOMVcombo_HV1.000_IHV1.000_CENV0.300\n')
          else:
            ofp.write('FASTP_FBASE       precomputeMWhiresCOMVcombo_HV1.000_IHV1.000_CENV0.000\n')
          if runslow == 1:
            ofp.write('USE_FASTP_COUNTS      0\n')
          else:
            ofp.write('USE_FASTP_COUNTS    2\n')
        if(use5002or5003 == 0):
          ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5002.txt\n')
          ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1_splits5_6.syswtheory\n')
          ofp.write('XIELL_REBINFNAME    bin1.txt.reformat.new\n')
          ofp.write('ELLMAX              2\n')
        else:
          ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5003.txt\n')
          ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1_splits5_6.syswtheory.dummy\n')
          ofp.write('XIELL_REBINFNAME    bin1.txt.reformat.new\n')
          ofp.write('ELLMAX              2\n')
      else:
        ## temporary!
        ## needed before had the full grid.
        #assert np.fabs(hvscale - 1.) < 2.0e-6
        #assert np.fabs(ihvscale - 1.) < 2.0e-6
        #assert np.fabs(cenvfrac - 0.) < 2.0e-6
        #ofp.write('FASTP_FBASE         precomputeMWhiresmaskedcombo_HV1.000_IHV1.000_CENV0.000\n')
        #ofp.write('USE_FASTP_COUNTS    2\n')
    
        if COMVopt == 0:
          if (np.fabs(hvscale - 1.0) < 1e-4 and \
             np.fabs(ihvscale - 1.0) < 1e-4 and \
             np.fabs(cenvfrac) < 1e-4) or (runslow == 1):
             ofp.write('FASTP_FBASE    precomputeMWhiresmaskedcombo_HV1.000_IHV1.000_CENV0.000\n')
             if runslow == 1:
               ofp.write('USE_FASTP_COUNTS      0\n')
             else:
               ofp.write('USE_FASTP_COUNTS    2\n')
          elif useprecp2 is not None:
            ofp.write('FASTP_FBASE        %s\n' % (useprecp2))
            ofp.write('USE_FASTP_COUNTS    2\n')
          else:  ## use the grid.
            ofp.write('FASTP_FBASE         gridv0masked/precomputeMWhiresmaskedcombo.precptable\n')
            assert runslow == 0
  #          if runslow == 1:
  #            ofp.write('USE_FASTP_COUNTS      0\n')
  #          else:
            ofp.write('USE_FASTP_COUNTS    3\n')
  
        else:  ##COMVopt = 1
          #### add this option later?
          #elif useprecp2 is not None:
          #  ofp.write('FASTP_FBASE        %s\n' % (useprecp2))
          #  ofp.write('USE_FASTP_COUNTS    2\n')
          assert np.fabs(cenvfrac) < 0.001 or np.fabs(cenvfrac - 0.3) < 0.001 or runslow == 1
          if np.fabs(cenvfrac - 0.3) < 0.001:
            ofp.write('FASTP_FBASE       precomputeMWhiresCOMVmaskedcombo_HV1.000_IHV1.000_CENV0.300\n')
          else:
            ofp.write('FASTP_FBASE       precomputeMWhiresCOMVmaskedcombo_HV1.000_IHV1.000_CENV0.000\n')
          if runslow == 1:
            ofp.write('USE_FASTP_COUNTS      0\n')
          else:
            ofp.write('USE_FASTP_COUNTS    2\n')
  
        if use5002or5003 == 0:
          ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5002bin2.txt.dummyrow\n')
          ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1fineMU_splits1_1.syswtheory.dummyrow\n')
          ofp.write('XIELL_REBINFNAME    bin1fineMU.txt.reformat.new\n')
          ofp.write('ELLMAX              2\n')
        else:
          ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5003bin2.txt.dummy\n')
          ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1fineMU_splits1_1.syswtheory.dummy\n')
          ofp.write('XIELL_REBINFNAME    bin1fineMU.txt.reformat.new\n')
          ofp.write('ELLMAX              2\n')
## end whichbox == 1 for xiell.

    if whichbox == 2:
      assert use5002or5003 == 1
      if maskedopt == 1:
        ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5003bin2.txt.dummy\n')
        ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1fineMU_splits1_1.syswtheory.dummy\n')
        ofp.write('XIELL_REBINFNAME    bin1fineMU.txt.reformat.new\n')
        ofp.write('ELLMAX              2\n')
  #      ofp.write('FASTP_FBASE           precomputePB00maskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n')
  ## moved up.
  #      ofp.write('FASTP_FBASE_WP        precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n' % (whichPB))
        if useprecp2 is None:
          if runslow == 1:
            ofp.write('USE_FASTP_COUNTS    0\n')
            ofp.write('FASTP_FBASE        precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n' % (whichPB))
          else:
            if np.fabs(ihvscale - 1.) < 0.0001:
              ofp.write('FASTP_FBASE       PB%02d/precomputePB%02dmaskedcombo_0.6452.precptable\n' % (whichPB,whichPB))
            elif np.fabs(ihvscale - 1.2) < 0.0001:
              ofp.write('FASTP_FBASE       PB%02dIHV1p2/precomputePB%02dmaskedcombo_0.6452.precptable\n' % (whichPB,whichPB))
            else:
              print 'unsupported IHV',ihvscale
              return None
            ofp.write('USE_FASTP_COUNTS    3\n')
        else:
          ofp.write('USE_FASTP_COUNTS    2\n')
          mystr = 'PB%02d' % whichPB
  ## make sure PB set in agreement with wp through whichPB!
          assert re.search(mystr,useprecp2)
          ofp.write('FASTP_FBASE        %s\n' % (useprecp2))
      else:  #maskedopt == 0
        assert runslow == 1
        ofp.write('USE_FASTP_COUNTS    0\n')
        ofp.write('FASTP_FBASE        precomputePB%02dmaskedcombo_0.6452_HV1.000_IHV1.000_CENV0.000\n' % (whichPB))

        if(use5002or5003 == 0):
          ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5002.txt\n')
          ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1_splits5_6.syswtheory\n')
          ofp.write('XIELL_REBINFNAME    bin1.txt.reformat.new\n')
          ofp.write('ELLMAX              2\n')
        else:
          ofp.write('XIELL_DATAFNAME     xi02NSdebiasedboss5003.txt\n')
          ofp.write('XIELL_ICOVFNAME     icovtotv7corr_b5000000_N200_rebin-bin1_splits5_6.syswtheory.dummy\n')
          ofp.write('XIELL_REBINFNAME    bin1.txt.reformat.new\n')
          ofp.write('ELLMAX              2\n')


  if whichbox == 0:
    assert(COMVopt == 0) ## otherwise, we just need to hard code paths to the COMV catalogs.
    Lbox = 2750.0
    ofp.write('HaloFileName     /home/howdiedoo/SOforL0work/packSOforL0openmpnewvel/SOforL0.concat\n')
    ofp.write('HaloDmFileName   /home/howdiedoo/SOforL0work/packSOforL0openmpnewvel/SOforL0.concat.halomembers\n')
    ofp.write('FOFHaloFileName  /home/howdiedoo/SOforL0work/packSOforL0openmp/SOconvertboundary/FOFnomatchfile.dat.cut\n')
    ofp.write('DmFileName       /home/howdiedoo/SOforL0work/packSOforL0openmp/L00_0.6452.dm\n')
    ofp.write('MASSFXNFNAME     massfxndlg10m0p01.dat\n')
  elif whichbox == 1:
    Lbox = 677.7
    ofp.write('HaloFileName     /home/howdiedoo/SOmaster/MWhiresplanck.halos\n')
    ofp.write('HaloDmFileName   /home/howdiedoo/SOmaster/MWhiresplanck.halomembers\n')
    ofp.write('FOFHaloFileName  doesnotexist\n')
    ofp.write('DmFileName       /home/howdiedoo/SOmaster/MWhires.dm\n')
    ofp.write('MASSFXNFNAME     massfxndlg10m0p01.dat.blerg.gack\n')
  elif whichbox == 2:
    Lbox = 1380.0
    ofp.write('HaloFileName     /home/howdiedoo/SOmaster/PB%02d_0.6452.halos\n' % (whichPB))
    ofp.write('HaloDmFileName   /home/howdiedoo/SOmaster/PB%02d_0.6452.halomembers\n' % (whichPB))
    ofp.write('FOFHaloFileName  doesnotexist\n')
    ofp.write('DmFileName       doesnotexist\n')
    ofp.write('MASSFXNFNAME     doesnotexist\n')

  else:
    print 'whichbox = ',whichbox,'not supported.'
    return None

  if not ifp:
    print 'set up %s file first, then rerun' % (parambase)
    return None
  for line in ifp:
    ofp.write('%s' % (line))
  ifp.close()

  ofp.write('HVSCALE    %f\n' % (hvscale))
  ofp.write('IHVSCALE    %f\n' % (ihvscale))
  ofp.write('CENVFRAC    %f\n' % (cenvfrac))
  ofp.write('BOX_SIZE    %f\n' % (Lbox))
  ofp.write('CENOPT     %d\n' % (cenopt))
  ofp.write('SATOPT     %d\n' % (satopt))
  ofp.write('FORCE_CENFORSAT   %d\n' % (force_cenforsat))
  ofp.write('PARAMFNAME   %s\n' % (hodnopath))
  ofp.write('OUTFILETAG   %s\n' % (outfiletag))
  ofp.write('FIT_LENSING   %d\n' % (whichfit[0]))
  ofp.write('FIT_MULTIPOLES   %d\n' % (whichfit[1]))
  ofp.write('FIT_WP   %d\n' % (whichfit[2]))
  ofp.close()

  mycmd = './runall %s > runalltmpout' % (batfname)
  print 'running runall:',mycmd
  os.system(mycmd)
  print 'run finished.  moving back to original directory.'
  os.chdir(mycwd)
  return outfiletag ## know where to get the result for plotting, etc after run.

def runchainmodel(celt,ftag,pbase,maskedopt=0,whichfit=[0,1,0],whichbox=1,cenopt=2,satopt=2, force_cenforsat=1,setihv0=0, setcenv0=0, setnosats=0, hvscalenew=None, ihvscalenew=None, cenvnew=None, runslow=0, COMVopt=0,use5002or5003=0,useprecp2=None,whichPB=0,writecat=None, aperp = 1.0, apar = 1.0):
  """
  If you want to turn off velocities for comparison, setihv0=1 will set ihvscale=0.
  setcenv0 will set cenvfrac=0
  setnosats will set Mcut = 2.e16.
  useprecp2 = None if no precp2 available, useprecp2 = FASTP_BASE of precp2 case.
  """
  logopt = 0
  if celt['M_min'] < 1e10:
    logopt = 1 

  if logopt == 1:
    Mmin = 10**(celt['M_min'])
    Mcut = 10**(celt['M_cut'])
    M1 = 10**(celt['M1'])
  else:
    Mmin = (celt['M_min'])
    Mcut = (celt['M_cut'])
    M1 = (celt['M1'])

  myihv = celt['ihvscale']
  mycenv = celt['cenvfrac']
  myhv = celt['hvscale']
  if setnosats == 1:
    Mcut = 2.e16
  if setihv0 == 1:
    myihv = 0.

  if setcenv0 == 1:
    mycenv = 0.

  if hvscalenew is not None:
    myhv = hvscalenew

  if ihvscalenew is not None:
    myihv = ihvscalenew

  if cenvnew is not None:
    mycenv = cenvnew


  outfiletag = runmodel(Mmin = Mmin, M1 = M1, alpha = celt['alpha'], \
           Mcut = Mcut, sigmalogM = celt['sigma_logM'],\
           hvscale = myhv, ihvscale = myihv, cenvfrac = mycenv,\
           whichbox=whichbox,cenopt=cenopt, satopt=satopt, force_cenforsat=force_cenforsat,\
           ftag=ftag,pbase=pbase,maskedopt=maskedopt,whichfit=whichfit, runslow=runslow, COMVopt = COMVopt,\
           use5002or5003=use5002or5003,useprecp2=useprecp2,whichPB=whichPB,writecat=writecat,aperp=aperp,apar=apar)

  if outfiletag is None:  return None
  return pbase + 'fits2data/'+outfiletag 

def subsamplechainslow(chainfname,colfname,nsubsample,pbase,use5002or5003=1):

  ccx = mcmcutils.chain(chainfname,colfname)
  wgtsum = int(ccx.chain['weight'][:].sum())
  cumwgt = ccx.chain['weight'][:].cumsum()
  xx = np.random.randint(0,wgtsum,nsubsample)
  xxi = np.zeros(len(xx))
  for i in range(len(xxi)):
    xxi[i] = np.where(cumwgt < xx[i])[0][-1]
    assert cumwgt[xxi[i]] < xx[i]
    if xxi[i] < nsubsample - 1:
      assert cumwgt[xxi[i]+1] >= xx[i]
  assert ((xxi >= 0) & (xxi < len(ccx.chain['weight'][:]))).all()
  print 'generated this many subsamples',len(xxi)
  for i in range(len(xxi)):
    ftag = chainfname.split('/chains/')[1].split('.chain')[0] + '_precp0_%06d_' % xxi[i]
    runchainmodel(ccx.chain[xxi[i]], pbase=pbase,maskedopt=1,whichfit=[0,1,1],runslow=1,ftag=ftag,use5002or5003=use5002or5003)
  #return xxi


if __name__ == '__main__':
  pass



