#!/usr/bin/env python

import subprocess
import sys
import os.path as osp

if (len(sys.argv)!=5):
    sys.exit("ABORT: filen, # of variations, # of configurations, force submit (0:false) needed")

force = 0

filen = sys.argv[1]
nr    = int(sys.argv[2])
nconf = int(sys.argv[3])
force = int(sys.argv[4])

for i in range(nr):
    for j in range(nconf):
        pass
        if( osp.isfile( "Smat."+filen+".%i.%i.dat"%(i,j) ) == False or force!=0):
            subprocess.call(["qsub", "submit.%i.%i.sh"%(i,j)])


