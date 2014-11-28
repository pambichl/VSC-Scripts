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
out_str = [open("qpeak.%i.xml"%i, "w") for i in range(len(pos_range))]
for line in in_str:
    line = line.strip() #entfernt alle trailing characters vorne und hinten
    for i in range(len(pos_range)):
        if line.startswith('<input> Abmessung.xml </input>'):
            line_out = line + '\n<param name="file">' + filen + '.qpeaks.%i</param>'%i
        else:
            line_out = line
        out_str[i].write(line_out)
        out_str[i].write("\n")
for s in np.ravel(out_str):
    s.close()

#Submit-Skripts erstellen
in_str = open("peaksubmit.sh", "r").readlines()
out_str = [open("qpeaksubmit.%i.sh"%i, "w") for i in range(len(pos_range))]
for line in in_str:
    line = line.strip()
    for i in range(len(pos_range)):
        if 'mpirun' in line:
            line = 'time mpirun -machinefile $TMPDIR/machines -np $NSLOTS solve_xml_mumps -i qpeak.%i.xml'%i
            out_str[i].write(line)
        else:
            out_str[i].write(line)
        out_str[i].write("\n")
for s in np.ravel(out_str):
    s.close()
                
#jobs ausfuehren         
for i in range(len(pos_range)):
    pass
    subprocess.call(["qsub.py", "qpeaksubmit.%i.sh"%i])
