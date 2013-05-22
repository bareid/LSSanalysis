#import numpy as N
import os
import re
import sys

def riemannheader(ofp,jobname,jobtime='12:00:00',ppn=1,fastqopt=0):

  ofp.write('#!/bin/bash\n')
  ofp.write('#PBS -N %s\n' % jobname)
  ofp.write('#PBS -l nodes=1:ppn=%d,walltime=%s\n' % (ppn,jobtime))
  ofp.write('#PBS -o %s.$PBS_JOBID.out\n' % (jobname))
  ofp.write('#PBS -e %s.$PBS_JOBID.err\n' % (jobname))
  if(fastqopt == 1):
    ofp.write('#PBS -q fast\n')
  ofp.write('#PBS -V\n')
  ofp.write('#\n')
  ofp.write('cd $PBS_O_WORKDIR\n')

def nerscheader(ofp,jobname,jobtime='23:59:59',mppwidth=24,qname='thruput'):
  ofp.write('#PBS -q %w\n' % (qname))
  ofp.write('#PBS -l mppwidth=%d\n' % (mppwidth))
  ofp.write('#PBS -l walltime=%s\n' % (jobtime))
  ofp.write('#PBS -j eo\n')
  ofp.write('#PBS -V\n')
  ofp.write('#\n')
  ofp.write('cd $PBS_O_WORKDIR\n')
  ofp.write('#\n')
  


if __name__ == '__main__':

  if(len(sys.argv) != 3):
    print 'Usage: makexijobscripts.py task whichmachine'
    print 'note this will break if you dont run it on whichmachine'
    print 'whichmachine = 0 [howdie]'
    print 'whichmachine = 1 [riemann]'
    print 'whichmachine = 2 [NERSC]'
    print 'edit file directly to change riemann/NERSC queue'
    print 'task = 0: DD, DR, RR counts needed for angular weight calculation'
    print 'task = 1: convert D/R counts to angweight and copy it over to other directories that need it'
    print 'task = 2: do xiell counts.'
    sys.exit(1)

  fastqopt = 0

  task = int(sys.argv[1])
  whichmachine = int(sys.argv[2])

  ## data dir locations.
  ## HH = howdie, RR = riemann, NN = nersc.

  dHH = "/home/howdiedoo/boss/tiledmockboss5002redo/"
  dRR = "/data/bareid/tiledmockboss5002redo/"
  dNN = "/scratch/scratchdirs/bareid/tiledmockboss5002redo/"

  ## working directory base.
  wHH = "/home/howdiedoo/boss/"
  wRR = "/home/bareid/boss/tiledmockboss5002redo/"
  wNN = "/global/u1/b/bareid/boss/"      

  ## nersc
  mppwidth = 24

  ## machine specific pre/post commands to the main C call.
  preambleHH = "nohup "
  preambleRR = ""
  preambleNN = "aprun -n 1 -N 1 -d %d " % mppwidth

  postambleHH = " & \n"
  postambleRR = "\n"
  postambleNN = "\n"

  dlist = [dHH, dRR, dNN]
  wdirlist = [wHH, wRR, wNN]
  prelist = [preambleHH, preambleRR, preambleNN]
  postlist = [postambleHH, postambleRR, postambleNN]

  ## data, random, fout, angweight pairs.
  cat012 = {'dataf':'cmass-boss5002sector-FBBRv2icoll012.dat','ranf':'randoms-boss5002-icoll012-vetoed.dat','outf':'outputtiledmockboss5002redo/cmass-boss5002sector-FBBRv2icoll012','angweight':None}
  catNN = {'dataf':'cmass-boss5002sector-FBBRv2-NN.dat','ranf':'randoms-boss5002-NN-vetoed.dat','outf':'outputtiledmockboss5002redo/cmass-boss5002sector-FBBRv2-NN','angweight':None}
  ## these are fed to angular clustering calculator only.
  catangtest = {'dataf':'cmass-boss5002sector-FBBRv2-ang.dat','ranf':'randoms-boss5002-ang-vetoed.dat','outf':'outputtiledmockboss5002redo/cmass-boss5002sector-FBBRv2-ang','angweight':'angweightboss5002redo.1pw'}
  catangzcut = {'dataf':'cmass-boss5002sector-FBBRv2-ang.dat','ranf':'randoms-boss5002-ang-vetoed.dat','outf':'outputtiledmockboss5002redo/cmass-boss5002sector-FBBRv2-angzcut','angweight':'angzcutweightboss5002redo.1pw'}

  ## these are fed upwards -- has the correct angular clustering upweighting file.
  catangfinal = {'dataf':'cmass-boss5002sector-FBBRv2-ang.dat','ranf':'randoms-boss5002-ang-vetoed.dat','outf':'outputtiledmockboss5002redo/cmass-boss5002sector-FBBRv2-ang','angweight':'angzcutweightboss5002redo.1pw'}

### need these for data, not the mocks though
#  catanghigh = {'dataf':,'ranf':,'outf':,'angweight':''}
#  catanglow = {'dataf':,'ranf':,'outf':,'angweight':''}

  ## jobs to do.
  ## job name, directory name, whichcat, fmtstring, nproc, jobtime
  ##

  ## depending on the inputs to the correlation fxn code, the string formats are slightly different.
  ## that's encoded in fmttype variable.  0 is std, 1 adds on the angular weight as a last argument.
  
  ## ang jobs for mocks.
  jobang1 = {'name': 'ang1', 'dname': 'zdistvXlogbinsompcleverLSangfaster', 'whichcat': catangtest, \
         'fmtstring': '%s./xiLSlogbins %s %s 4 4 %d %s > %s%s',  'fmttype':0, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  jobang1zcut = {'name': 'ang1zcut', 'dname': 'zdistvXlogbinsompcleverLSangfaster', 'whichcat': catangzcut, \
         'fmtstring': '%s./xiLSlogbins %s %s 4 4 %d %s 0.43 0.7 0 > %s%s ', 'fmttype':0, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  ## for these two, input list is (preamble, dfilename, outfname, DRval, rfilename, txtoutfname, postamble)
  ## xiell jobs.
  jobxiell012 = {'name': 'xiell012_boss5002', 'dname': 'zdistvXlogbinsompcleverLSsmallscale/', 'whichcat': cat012, \
         'fmtstring': '%s./xiLSmultipolesnoweight %s %s 2 %d %s > %s%s ', 'fmttype':0, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  jobxiellNN = {'name': 'xiellNN_boss5002', 'dname': 'zdistvXlogbinsompcleverLSsmallscale/', 'whichcat': catNN, \
         'fmtstring': '%s./xiLSmultipolesnoweight %s %s 2 %d %s > %s%s ', 'fmttype':0, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  jobxiellang = {'name': 'xiellang_boss5002', 'dname': 'zdistvXlogbinsompcleverLSsmallscale/', 'whichcat': catangfinal, \
         'fmtstring': '%s./xiLSmultipolesangweight %s %s 2 %d %s %s > %s%s ', 'fmttype':1, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}

  ## wp jobs.
  jobwp012 = {'name': 'wp012_boss5002', 'dname': 'zdistvXlogbinsompcleverLSsmallscale/', 'whichcat': cat012, \
         'fmtstring': '%s./xiLSwpnoweight %s %s 2 %d %s > %s%s ', 'fmttype':0, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  jobwpNN = {'name': 'wpNN_boss5002', 'dname': 'zdistvXlogbinsompcleverLSsmallscale/', 'whichcat': catNN, \
         'fmtstring': '%s./xiLSwpnoweight %s %s 2 %d %s > %s%s ', 'fmttype':0, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  jobwpang = {'name': 'wpang_boss5002', 'dname': 'zdistvXlogbinsompcleverLSsmallscale/', 'whichcat': catangfinal, \
         'fmtstring': '%s./xiLSwpangweight %s %s 2 %d %s %s > %s%s ', 'fmttype':1, 'DRlist':[1,2,3], \
         'nprocRR': 8, 'nprocNN':24, 'jobtimeRR':'23:59:00','jobtimeNN':'23:59:00'}
  
         

  ## task0: run code to get w(theta), from which to derive angular weights.
  ## hard code times and np here.

  ddir, workingdir, preamble, postamble = dlist[whichmachine], wdirlist[whichmachine], prelist[whichmachine], postlist[whichmachine]
  if(task==0):
    
    ## jobs in this list must have the same format.
    joblist = [jobang1, jobang1zcut]

  if(task == 1):
    #print 'task = 1: convert D/R counts to angweight and copy it over to other directories that need it'
    print 'to be written!'
    joblist = []
    sys.exit(1)

  ## task 2: run xiell.
  if(task == 2):
    joblist = [jobxiell012, jobxiellNN, jobxiellang]

  ## task 3: run wp
  if(task == 3):
    joblist = [jobwp012, jobwpNN, jobwpang]

  ## task 4:


  for jj in joblist:
    filesh = workingdir+jj['dname']+'/'+jj['name']+'.sh'
#      print 'sh is going here:',filesh
    ofpsh = open(filesh,'w')
    dfilename = ddir+jj['whichcat']['dataf']
    rfilename = ddir+jj['whichcat']['ranf']
    outfilename = workingdir+jj['dname']+'/'+jj['whichcat']['outf']
    if(whichmachine == 0):
      ofpsh.write('#!/bin/bash\n')
    if(whichmachine == 1):
      riemannheader(ofpsh,jj['name'],jobtime=jj[jobtimeRR], ppn=jj['nprocRR'], fastqopt=fastqopt)
    if(whichmachine == 2):
      nerscheader(ofpsh,jj['name'],jobtime=jj[jobtimeNN], mppwidth=jj['nprocNN'], qname='thruput')
    for DR in jj['DRlist']:
      txtout = jj['name']+'_'+str(DR)
      if(jj['fmttype'] == 0):
        mystr = jj['fmtstring'] % (preamble, dfilename,outfilename,DR,rfilename,txtout,postamble)
      if(jj['fmttype'] == 1):
        mystr = jj['fmtstring'] % (preamble, dfilename,outfilename,DR,rfilename, jj['whichcat']['angweight'], txtout,postamble)
      ofpsh.write( '%s' % (mystr)) 
    ofpsh.close()
    mystr = 'chmod +x %s' % (filesh)
    os.system(mystr)

    

