#!/usr/bin/env python

import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

if (len(sys.argv)!=5):
    sys.exit("ABORT: filen, # of calculations, # of slices and # parameter variations needed")

filen = sys.argv[1]
nr = int(sys.argv[2])
slicnum = int(sys.argv[3])
nparam = int(sys.argv[4])

ntimes = 12
freqs = range(0,nparam)

vals = []
for i in range (0,nr):
    calcs = []
    for j in range (0,slicnum):
        delaytimes = np.loadtxt("avedt."+filen+".%i.%i.dat"%(i+1,j+1))
        calcs.append(delaytimes)
    vals.append(calcs)

for i in range (0,nr):
    timeskinds  = []
    meankinds   = []
    stddevkinds = []
    for j in range (0,ntimes):
        timesscan   = []
        meanscan    = []
        stddevscan  = []
        for k in range (0,nparam):
            timesslices = []
            for l in range (0,slicnum):
                timesslices.append(vals[i][l][k + j*nparam][1])
            timesscan.append(timesslices)

            meanscan.append(np.mean(timesslices))
            stddevscan.append(np.std(timesslices))

        timeskinds.append(timesscan)
        meankinds.append(meanscan)
        stddevkinds.append(stddevscan)


    plt.figure(1)

    plt.errorbar(freqs, meankinds[0], stddevkinds[0], fmt='D-r', label=r'$\tau_{tot}$')

    #polycoeffs = np.polyfit(freqs, meankinds[0], 5)
    #poly = np.poly1d(polycoeffs)
    #plt.plot(freqs, poly(freqs), '--y', label=r'$\tau_{tot} \, \rm{fit}$', lw=2.0)

    Weyl = [10*2*10 * 0.96 / (2*Pi)] * nparam
    plt.plot(freqs, Weyl, '-g', label=r"$\rm{Weyl-estimate}$", lw=2.0)

    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'delay times')
    plt.axis([0,50, 0.,50.])
    # inches tight funkt anscheinend am vsc nicht!
    plt.savefig("NQtimes."+filen+".%i.png"%i, bbox_inches='tight')
    plt.clf()

    plt.figure(2)

    plt.errorbar(freqs, meankinds[1], stddevkinds[1], fmt='d-r', label=r'$\tau_{11}$')
    plt.errorbar(freqs, meankinds[2], stddevkinds[2], fmt='d-g', label=r'$\tau_{t}$')
    plt.errorbar(freqs, meankinds[3], stddevkinds[3], fmt='d-b', label=r'$\tau_{r}$')
    polycoeffs = np.polyfit(freqs, meankinds[1], 5)
    poly = np.poly1d(polycoeffs)
    plt.plot(freqs, poly(freqs), '--y', label=r'$\tau_{11} \, \rm{fit}$', lw=2.0)
    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'delay times')
    plt.savefig("NLefttimes."+filen+".%i.png"%i, bbox_inches='tight')
    plt.clf()

    plt.figure(3)

    plt.errorbar(freqs, meankinds[-1], stddevkinds[-1], fmt='d-r', label=r'$T$')
    plt.errorbar(freqs, meankinds[-2], stddevkinds[-2], fmt='d-g', label=r'$R$')
    plt.plot(freqs, np.array(meankinds[-1]) + np.array(meankinds[-2]), 'd--y', label=r'$T+R$')
    #polycoeffs = np.polyfit(freqs, meankinds[1], 5)
    #poly = np.poly1d(polycoeffs)
    #plt.plot(freqs, poly(freqs), '--y', label=r'$\tau_{11} \, \rm{fit}$', lw=2.0)
    plt.axis((freqs[0],freqs[-1], 0., 1.2))
    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'$\rm{Transmittivity/Reflectivity}$')
    plt.savefig("T+R."+filen+".%i.png"%i, bbox_inches='tight')
    plt.clf()


sp.call(["sh", "calcav.sh"])

    

        



  

