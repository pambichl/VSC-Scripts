#!/usr/bin/env python

import subprocess
import sys
import numpy as np
import os
import random

if (len(sys.argv)!=6):
    sys.exit("ABORT: parameters \nfilen, \n# of freq-variations, \n# of configurations, \nstart wavenumber [Pi/W], \nend wavenumber [Pi/W] \nneeded")

reffilen = str(sys.argv[1])   # namestring of calculation
nfreq    = int(sys.argv[2])   # number of radius variations
nconf    = int(sys.argv[3])   # number of configurations to be averaged over
startf   = float(sys.argv[4]) # start k in units of Pi/system width
endf     = float(sys.argv[5]) # end k in units of Pi/system width

#Vorlageninput zeilenweise in ein array einlesen
in_str  = open("input.xml", "r").readlines()
out_str = [[open("input.%i.%i.xml" % (i,j), "w") for j in range(nconf)] for i in range(nfreq)]

freqs = np.linspace(startf,endf,nfreq)

#Vorlageninput zeilenweise durchgehen und entsprechende Zeilen ersetzen
for line in in_str:

    line = line.strip() #entfernt alle trailing characters vorne und hinten

    for i in range(nfreq):
        for j in range(nconf):

            if line.startswith('<param name="modes">'):
                line_out = '<param name="modes">%.5f</param>' % freqs[i]
                out_str[i][j].write(line_out)
            elif line.startswith('<param name="file">av'):
                line_out = '<param name="file">' + reffilen + '.%i.%i</param>' % (i,j)
                out_str[i][j].write(line_out)
            elif line.startswith('<param name="srand">'):
                line_out = '<param name="srand">%i</param>' % random.randint(1,10000000)
                out_str[i][j].write(line_out)
            elif line.startswith('<param name="output">obstacles.dat'):
                line_out = '<param name="output">obstacles.' + reffilen + '.%i.%i.dat</param>' % (i,j)
                out_str[i][j].write(line_out)
            elif line.startswith('<param name="file">obstacles.dat'):
                line_out = '<param name="file">obstacles.' + reffilen + '.%i.%i.dat</param>' % (i,j)
                out_str[i][j].write(line_out)
            else:
                line_out = line
                out_str[i][j].write(line_out)

            out_str[i][j].write("\n")

for s in np.ravel(out_str):
    s.close()

#Submit-Skripts erstellen
in_str = open("submit.sh", "r").readlines()
out_str = [open("submit.0.%i.sh" % j, "w") for j in range(nconf)]

for line in in_str:
    line = line.strip()

    for j in range(nconf):

        if 'mpirun' in line:
            line = ''
            for i in range(nfreq):
                line += ('time mpirun -machinefile $TMPDIR/machines -np $NSLOTS solve_xml_mumps -i input.%i.%i.xml' % (i,j) + '\n' +
                         'rm pic.'+reffilen+'.*' + '\n')                
            out_str[j].write(line)
        else:
            out_str[j].write(line)

        out_str[j].write("\n")

for s in np.ravel(out_str):
    s.close()


