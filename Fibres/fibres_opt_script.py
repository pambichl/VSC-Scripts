#!/usr/bin/env python

import numpy as np

import subprocess
import os
from optparse import OptionParser


Pi = np.pi
I = 0. + 1.J


# read in time parameters from command line
parser = OptionParser()
parser.add_option("-f", "--filen",
                  action="store",
                  type = "str",
                  dest = "filen")
parser.add_option("-t", "--tries",
                  action="store",
                  type = "int",
                  dest = "tries")
options = parser.parse_args()[0]
filen = options.filen
tries = options.tries
if (not filen): print "ERROR: no filen given"
if (not tries): print "WARNING: no # of tries given, assuming 10"


dummy = open("dummy.txt", "w")
for i in range(tries): 
   print "run %4i of %4i" % (i+1,tries)

   randseed = np.random.randint(0,1000000)
   tag = "rand.auto.%04i" % (i+1)

   if (not os.path.exists(filen+".opt.log")):
      with open(filen+".opt.log", "w") as fstream:
         fstream.write("%4i %7i" % (i,randseed))
   else:
      with open(filen+".opt.log", "a") as fstream:
         fstream.write("%4i %7i" % (i,randseed))

   subprocess.call(["fibres_optimize.py",
                    "-n", "150000",
                    "-t", "2",
                    "-s", "0.5",
                    "-e", "%i" % randseed,
                    "--auto", "1"
                    ], stdout=dummy)
   subprocess.call(["fibres_opt_mov.py",
                    "--tmin", "-1.0",
                    "--tmax", "20.0",
                    "-n", "300",
                    "-t", tag,
                    "--auto", "1"
                    ], stdout=dummy)
   subprocess.call(["rm movie/*auto*.jpg"], shell=True)


dummy.close()
subprocess.call(["rm", "dummy.txt"])



