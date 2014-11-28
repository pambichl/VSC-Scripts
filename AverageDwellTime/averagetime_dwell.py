#!/usr/bin/env python

import numpy as np
import scipy as sp
import scipy.optimize as spopt
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

import warnings
import subprocess

from Scattering import scattering2D as scat
from Calc import fitfuncs as ff
from Calc import transmission as trans
from Calc import localization as local


if (len(sys.argv)!=6):
    sys.exit("ABORT: parameters filen, # of radius variations, # of configurations,  # of frequency variations, length of system [W] needed")

filen = str(sys.argv[1]) # namestring of calculation
nr    = int(sys.argv[2]) # number of radius variations
nconf = int(sys.argv[3]) # number of configurations to be averaged over
nfreq = int(sys.argv[4]) # number of frequency steps scanned over
L = float(sys.argv[5]) # length of system in units of width

data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/scatterdata/"

energs = scat.read_S(data_direc, filen+".0.0")[2]
if (len(energs)==0):
    kmin = 3.1 * np.pi / 10.
    kmax = 4.9 * np.pi / 10.
    dk = 0.5 * (kmin+kmax) * 10**(-7)
else:
    kvals = np.sqrt(2*energs)
    dk = (kvals[-1] - kvals[-3]) / 2.0

dims = scat.read_S(data_direc, filen+".0.0")[1]/2
nin_max = np.max(dims)

def read_single_S_list(filen, i, j):
    filen_i = filen+".%i.%i" % (i,j)
    print "reading file: "+filen_i
    try:
        data = scat.read_S(data_direc, filen_i)
        return data
    except Exception:
        if 'problems' in globals():
            problems.append(filen_i) 
        print "WARNING: problem reading "+filen_i
        return np.zeros((1, nfreq, 2*nin_max, 2*nin_max), dtype='complex')

problems = []
S = np.array([[read_single_S_list(filen, i ,j)[0] for j in range(nconf)] for i in range(nr)])
print 'WARNING: problems encountered in files', problems
for prob in problems:
    subprocess.call(["rm", data_direc+"Smat."+prob+".dat"])

ni = 0.0001
gamma = ni*2*energs[1]

q_mean  = np.array([[scat.calc_average_dwelltime(S[i,j],gamma) for j in range(nconf)] for i in range(nr)])
means   = np.array([[np.mean(q_mean[i,:,j]) for j in range(nfreq)] for i in range(nr)]).real
stddevs = np.array([[np.std(q_mean[i,:,j])  for j in range(nfreq)] for i in range(nr)]).real



### Plots ###


plot_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"

try:
    obstacles = np.loadtxt(plot_direc + 'potential.dat', usecols=(1,))
    obst_mask = abs(obstacles) == 0
    obst_ratio = float(np.shape(obstacles[obst_mask])[0]) / float(np.shape(obstacles)[0])
    print "filling factor of obstacles: %5f" % (1.0-obst_ratio)
except:
    obst_ratio = 0.94
    print "WARNING: automatically assuming filling fraction to be 0.06"

for mean, stddev, indx in zip(means, stddevs, range(nr)):

    plt.errorbar(range(nfreq), mean, stddev, fmt='D-r', label=r'$\tau_{tot}$')

    # Weyl estimate for dwell time: A / 2Pi = W * L * (1-fill) / 2Pi
    Weyl = [10*L*10 * obst_ratio / (2*Pi)]*nfreq
    plt.plot(range(nfreq), Weyl, '--g', label=r"$\rm{Weyl-estimate}$", lw=2.0)

    mean_mean = np.mean(mean)
    plt.plot(range(nfreq), [mean_mean]*nfreq, '--', color='orange', label=r"$\rm{mean}$", lw=2.0)
    print "Weyl estimate, mean-value:", Weyl[0], mean_mean

    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'$\left\langle \rm{Tr}\left( Q \right) \right\rangle$')
    plt.ylim(0., 2*Weyl[0])
    plt.savefig(plot_direc+"AverageTime."+filen+".%i.png" % indx, bbox_inches='tight')
    plt.clf()

    

        



  

