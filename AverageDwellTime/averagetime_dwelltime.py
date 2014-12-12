#!/usr/bin/env python

import numpy as np
import scipy as sp
import scipy.optimize as spopt
import scipy.integrate as spint
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

import warnings
import subprocess

from Scattering import scattering2D as scat
from Utils import utils as ut

from Calc import fitfuncs as ff
from Calc import transmission as trans
from Calc import localization as local

imag_i = 1.0J

if (len(sys.argv)!=4):
    sys.exit("ABORT: parameters filen, # of frequency steps, # of configurations needed.")

filen = str(sys.argv[1]) # namestring of calculation
nfreq = int(sys.argv[2]) # number of frequencies scanned
nconf = int(sys.argv[3]) # number of configurations to be averaged over

data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"
plot_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"

def read_single_Qd_matrix(filen, w, j):
    '''
    Read in single Q_d matrix from file.
    '''
    global problems
    
    filen_i = "Qdmat." + filen + ".%i.%i.dat" % (w,j)
    location = data_direc + "dwelltimedata/" + filen_i
    print "reading file: " + filen_i
    try:
        f = open(location, 'r')
        energ = f.readline()
        dim = int(f.readline())
        f.close()
        data = np.loadtxt(location, dtype='float', skiprows=2, unpack=True)
        Qd = data[2] + imag_i*data[3]
        Qd = Qd.reshape((dim,dim))
        return Qd
    except:
        if 'problems' in globals():
            problems.append(filen_i) 
        print "WARNING: problem reading "+filen_i
        return 0.0

problems = []
data = np.array([[read_single_Qd_matrix(filen, w, j) for j in range(nconf)] for w in range(nfreq)])
print 'WARNING: problems encountered in files', problems

modes_min = 10.1
modes_max = 14.9
lead_width = 10.
kvals = np.linspace(modes_min, modes_max, nfreq) * np.pi / lead_width

print np.shape(data)
single_mean_q = np.array([[np.mean(np.diag(data[w,j].real)) for j in range(nconf)] for w in range(nfreq)])
print single_mean_q[0,0]
print np.trace(data[0,0].real)/20.
mean_q_list = kvals * np.mean(single_mean_q, axis=1)
stddev = kvals * np.std(single_mean_q, axis=1)
mean_q = np.mean(mean_q_list)


#################
###   Plot   ####
#################

try:
    obstacles = np.loadtxt(plot_direc + 'potential.dat', usecols=(1,))
    obst_mask = abs(obstacles) == 0
    obst_ratio = float(np.shape(obstacles[obst_mask])[0]) / float(np.shape(obstacles)[0])
    print "filling factor of obstacles: %5f" % (1.0-obst_ratio)
except:
    obst_ratio = 0.94
    print "WARNING: automatically assuming filling fraction to be 0.06"

plt.errorbar(range(nfreq), mean_q_list, stddev, fmt='D-r', label=r'$\tau_{tot}$')

L = 2.5

# Weyl estimate for dwell time: A / 2Pi = W * L * (1-fill) / 2Pi
Weyl = [10*L*10 * obst_ratio / (2*Pi)]*nfreq
plt.plot(range(nfreq), Weyl, '--g', label=r"$\rm{Weyl-estimate}$", lw=2.0)

plt.plot(range(nfreq), [mean_q]*nfreq, '--', color='orange', label=r"$\rm{mean}$", lw=2.0)
print "Weyl estimate, mean-value:", Weyl[0], mean_q

plt.legend( bbox_to_anchor=(1.,1.) )
plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
plt.ylabel(r'$\left\langle \tau_W \right\rangle$')
plt.ylim(0., 2*Weyl[0])
plt.savefig(plot_direc+"AverageTime.DwellTime."+filen+".png", bbox_inches='tight')
plt.clf()

