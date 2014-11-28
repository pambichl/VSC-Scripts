#!/usr/bin/env python

import subprocess
import sys

if (len(sys.argv)!=3):
    sys.exit("ABORT: two jobids needed")

startjob = int(sys.argv[1])
endjob   = int(sys.argv[2])

for i in range(startjob,endjob+1):
    subprocess.call(["qdel", "%i" % i])
