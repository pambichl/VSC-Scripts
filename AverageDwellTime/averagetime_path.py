#!/usr/bin/env python

import numpy as np
import scipy.optimize as spopt
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

from Scattering import scattering2D as scat
from Utils import utils as ut

imag_i = 1.j


if (len(sys.argv)!=2):
    sys.exit("ABORT: filen-string needed")
filen = str(sys.argv[1]) # namestring of calculation
data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen

par = ut.read_input(data_direc)
try:
    nr = int(par['nr']) # number of radius variations
    nconf = int(par['nconf']) # number of configurations to be averaged over
    nfreq = int(par['nfreq']) # number of frequency steps scanned over
    lead_width = float(par['lead_width'])
    nyout = int(par['nyout'])
except KeyError:
    raw_input("WARNING: parameter missing in pinput.dat")

#energs = scat.read_S(data_direc, filen+".0.0")[2]
#kvals = np.sqrt(2*energs)
#dk = (kvals[-1] - kvals[-3]) / 2.0

dims = np.take(scat.read_S(data_direc + "/scatterdata/", filen + ".0.0")[1]/2, np.arange(1,3*nfreq,3))
nin_max = np.max(dims)

def read_state(filen, i, j, w, m):
    filen_i = filen + ".%i.%i.%04i.%04i" % (i,j,w,m)
    print "reading state: "+filen_i
    try:
        data = np.loadtxt(data_direc + "/states/" + "state." + filen_i + ".streu.dat",
                          dtype ='float', skiprows=1, usecols=(1,2), unpack=True)
        state = data[0] + imag_i * data[1]
        return state
    except Exception:
        print "WARNING: problem reading state "+filen_i
        return np.array(0.)

def read_refind(filen, i, j):
    filen_i = filen + ".%i.%i" % (i,j)
    print "reading refractive index: "+filen_i
    try:
        refind = np.loadtxt(data_direc + "/refindices/" + "obstacles." + filen_i + ".dat",
                          dtype ='float', usecols=(1,), unpack=True)
        mask = refind < 0.
        refind[mask] = 0.
        return np.sqrt(refind+1)
    except Exception:
        print "WARNING: problem reading refractive index "+filen_i
        return np.array(1.)

sys_size = read_refind(filen, 0, 0).size
tot_size = read_state(filen, 0,0,1,0).size
start_ind = (tot_size - sys_size) / 2
end_ind = start_ind + sys_size

time_ind = 2
dx = lead_width / (nyout+1)
times = np.zeros((nr,nconf,nfreq,nin_max), dtype='float')
for i in range(nr):
    for j in range(nconf):
        refind = read_refind(filen, i,j)
        for w in range(nfreq):
            for m in range(dims[w]):
                state = read_state(filen, i,j,3*w+1,m)[start_ind:end_ind]
                if (~np.iscomplex(state)).any(axis=0):
                    print "WARNING: nan encountered in state."+ filen + ".%i.%i.%04i.%04i" % (i,j,w,m)
                    state = state[np.iscomplex(state)]
                    refind = refind[np.iscomplex(state)]
                time = (np.vdot(state, refind**(time_ind) * state)).real * dx**2 / 2.0
                times[i,j,w,m] = time

qmean = np.sum(times, axis=-1)
means = np.array([[np.mean(qmean[i,:,w]/dims[w])
                   for w in range(nfreq)]
                  for i in range(nr)])
stddevs = np.array([[np.std(qmean[i,:,w]/dims[w])
                     for w in range(nfreq)]
                    for i in range(nr)])

plot_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"

for mean, stddev, indx in zip(means, stddevs, range(nr)):

    plt.errorbar(range(nfreq), mean, stddev, fmt='D-r', label=r'$\tau_{tot}$')
    
    # Weyl estimate for dwell time: A / 2Pi = W * L * (1-fill) / 2Pi
    Weyl = [10*2.5*10 * 0.90 / (2*Pi)]*nfreq
    plt.plot(range(nfreq), Weyl, '--g', label=r"$\rm{Weyl-estimate}$", lw=2.0)
    
    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'$\left\langle \rm{Tr}\left( Q \right) \right\rangle$')
    plt.savefig(plot_direc+"AveragePath.%i"%time_ind + "." + filen + ".%i.png"%indx, bbox_inches='tight')
    plt.clf()
 


  

