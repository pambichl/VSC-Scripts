#!/usr/bin/env python

import sys
import subprocess
import os

import numpy as np

from Utils import utils as ut

par = ut.read_input('./')
try:
   filen = par['filen'] 
   modes_min = float(par['modes_min'])
   modes_max = float(par['modes_max'])
   nin = int(par['refmodes']) # nin = refmodes
   ptc = int(par['ptc'])
except KeyError:
    raw_input("WARNING: parameter missing in pinput.dat")

nin_Max = int(0.5 * (modes_max + modes_min)) # nin_Max = n_open_modes

pos_range = (0, int(0.33*nin-0.5*ptc), int(0.5*nin-0.5*ptc))#, int(0.67*nin-0.5*ptc), nin-ptc)

in_str  = open("peak.xml", "r").readlines()
out_str = open("ypeak.xml", "w")
for line in in_str:
    line = line.strip()
    if line.startswith('<input> Abmessung.xml </input>'):
        line_out = line + '\n<param name="file">' + filen + '.ypeaks</param>'
    else:
        line_out = line
    out_str.write(line_out)
    out_str.write("\n")
out_str.close()

#Submit-Skripts erstellen
in_str = open("peaksubmit.sh", "r").readlines()
out_str = open("ypeaksubmit.sh", "w")
for line in in_str:
    line = line.strip()
    if 'mpirun' in line:
        line = 'time mpirun -machinefile $TMPDIR/machines -np $NSLOTS solve_xml_mumps -i ypeak.xml'
        out_str.write(line)
    else:
        out_str.write(line)
    out_str.write("\n")
out_str.close()

#job ausfuehren 
subprocess.call(["qsub.py", "ypeaksubmit.sh"])
