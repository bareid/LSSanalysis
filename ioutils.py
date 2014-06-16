import os
import re
import sys
import numpy as np
import struct

def readparamfile(fields,fname,sepchar='='):

  ifp = open(fname,'r')
  mydict = {}
  for f in fields:
    mydict[f] = None

  for line in ifp:
    elts = line.split(sepchar)
    ## remove blanks
    for jndx in range(elts.count('')):
      elts.remove('')
    assert len(elts) == 2
    myfield = elts[0].strip()
    if myfield in mydict:
      mydict[myfield] = elts[1].strip()

  return mydict

  ifp.close()

def compareascii(fname1,fname2,plist,comparetype=0,printopt=0):
  """
  comparetype specifies how to compare to files of numbers.
  plist contains the parameters of that comparison.
  if plist needs only one value, send a number works.
  comparetype = 0 checks if absolute values are less than plist[0]
  printopt = 1 will print first 100 discrepant lines.
  returns -1 if file structure does not match
  otherwise returns number of mismatched lines.
  """
  if type(plist) is not list:
    plist = [plist]

  if(comparetype != 0):
    print 'comparetype not defined, returning -1'
    return -1

  ifp1 = open(fname1,'r')
  ifp2 = open(fname2,'r')
  cnt = 0
  badmatch = 0
  for line1 in ifp1:
    line2 = ifp2.readline()
    if line2 == '': ## end of file!
      print 'files have different lengths!'
      return -1
    cnt += 1
    a1 = np.array(line1.split(),dtype='float')
    a2 = np.array(line2.split(),dtype='float')
    if len(a1) != len(a2):
      print 'diff number of columns',len(a1),len(a2),cnt
      return -1
    if(comparetype == 0):
      if(np.fabs(a1 - a2) > plist[0]).any():
        badmatch += 1
        if(printopt == 1 and badmatch <= 100):  
          print line, line2

  ifp1.close()
  ifp2.close()
  return badmatch

def comparebinary(fname1,fname2,plist,formatstring,bytesperstruct,comparetype=0,printopt=0):
  """
  ** ASSUMES number of lines is given first as an int. **
  ** FOLLOWED by a list of structs unpacked with formatstring **
  formatstring for subsample files is 'ffffffi', which is 28 bytes.
  comparetype specifies how to compare to files of numbers.
  plist contains the parameters of that comparison.
  if plist needs only one value, send a number works.
  comparetype = 0 checks if absolute values are less than plist[0]
  printopt = 1 will print first 100 discrepant lines.
  returns -1 if file structure does not match
  otherwise returns number of mismatched lines.
  """
  if type(plist) is not list:
    plist = [plist]

  if(comparetype != 0):
    print 'comparetype not defined, returning -1'
    return -1

  ifp1 = open(fname1,'rb')
  ifp2 = open(fname2,'rb')

  ## nevermind, not in my std header.
#  x1 = ifp1.read(4)
#  N1 = (struct.unpack("i",x1))[0]
#  x2 = ifp2.read(4)
#  N2 = (struct.unpack("i",x2))[0]
#  if(N1 != N2):
#    print 'files have different lengths!',N1,N2
#    return -1
  badmatch = 0

#  for i in range(N1):  ## read subsample particles one at a time.
  while(0==0):
    x1 = ifp1.read(bytesperstruct)
    x2 = ifp2.read(bytesperstruct)
    if not x1:
      break
    if not x2:
      break
    a1 = np.array(struct.unpack(formatstring,x1))
    a2 = np.array(struct.unpack(formatstring,x2))
    if(comparetype == 0):
      if(np.fabs(a1 - a2) > plist[0]).any():
        badmatch += 1
        if(printopt == 1 and badmatch <= 100):  
          print a1
          print a2
          print '    '

  ifp1.close()
  ifp2.close()
  return batmatch


if __name__ == '__main__':

  myfields = ['particle_minimum','density_threshold','DMfilebase','density_file','Omega_m','density_min','search_radius','input_file_type','nfiles','dendir','outfilebase','Npcuberoot','Lbox']
  mydict = readparamfile(fields=myfields,fname="L0lowz.SOini")
  print mydict
  


