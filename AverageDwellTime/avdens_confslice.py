#!/usr/bin/env python

import subprocess
import numpy
import os
import random

#Name der Referenzrechnung
reffilen = "av12"
#Anzahl der Rechnungen
nr = 8
#Anzahl der Konfigurationsvariationen und
#in wieviele Intervalle diese unterteilt werden sollen
nconf = 100
slicnum = 100
if nconf%slicnum != 0:
    print ("WARNING: number of frequencies not divisible by number of slices")
#filling factors der Streuer (in % der Flaeche der Streuregion)
startfill = 0.01; endfill = 0.08

#Vorlageninput zeilenweise in ein array einlesen
instr  = open("input.xml", "r").readlines()

#leere dictionaries anlegen
outstr = {}
radii  = {}

#neue keys und Werte hinzufuegen
#outstr ist dictionary von lists
for i in range(1,nr+1):
    radii[i]  = numpy.linspace(startr,endr,nr)[i-1]
    outstrcalcs = {}
    for j in range(1,slicnum+1):
        outstrcalcs[j] = open("input.%i.%i.xml" % (i,j), "w")
    outstr[i] = outstrcalcs
 
#Vorlageninput zeilenweise durchgehen und entsprechende Zeilen ersetzen
for line in instr:

    line = line.strip() #entfernt alle trailing characters vorne und hinten

    if line.startswith('<param name="NR">'):
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                line = '<param name="NR">%i</param>' % int(nconf/slicnum)
                outstr[i][j].write(line)
    elif line.startswith('<param name="fillmin">'):
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                line = '<param name="fillmin">%.3f</param>' % fillfacs[i]
                outstr[i][j].write(line)
    elif line.startswith('<param name="fillmax">'):
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                line = '<param name="fillmax">%.3f</param>' % fillfacs[i]
                outstr[i][j].write(line)
    elif line.startswith('<param name="file">'):
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                line = '<param name="file">' + reffilen + '.%i.%i</param>' % (i,j)
                outstr[i][j].write(line)
    elif line.startswith('<param name="srand">'):
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                line = '<param name="srand">%i</param>' % random.randint(1,10000000)
                outstr[i][j].write(line)
    else:
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                outstr[i][j].write(line)
                
    for i in range(1,nr+1):
        for j in range(1,slicnum+1):
            outstr[i][j].write("\n")

for i in range(1,nr+1):
    for j in range(1,slicnum+1):
        outstr[i][j].close()

#Submit-Skripts erstellen
instr = open("submit72.sh", "r").readlines()

for i in range(1,nr+1):
    for j in range(1,slicnum+1):
        outstr[i][j] = open("submit72.%i.%i.sh" % (i,j), "w")

for line in instr:
    line = line.strip()

    if line == 'time mpirun -q 36000 -m $TMPDIR/machines -np $NSLOTS ./solve_xml_mumps':
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                line = 'time mpirun -q 36000 -m $TMPDIR/machines -np $NSLOTS ./solve_xml_mumps -i input.%i.%i.xml' % (i,j)
                outstr[i][j].write(line)
    else:
        for i in range(1,nr+1):
            for j in range(1,slicnum+1):
                outstr[i][j].write(line)

    for i in range(1,nr+1):
        for j in range(1,slicnum+1):
            outstr[i][j].write("\n")
                

#jobs ausfuehren         
#for i in range(1,nr+1):
#    for j in range(1,slicnum+1):
#        outstr[i][j].close()
#        subprocess.call(["qsub.py", "submit72.%i.%i.sh" % (i,j)])

