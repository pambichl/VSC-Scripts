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

imag_i = 1.0J

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

##### for av42
#k = 10.5 * np.pi / 10.
#dk = 0.5 * k * 10.**(-10)
####

dims = np.take(scat.read_S(data_direc, filen+".0.0")[1]/2, np.arange(1,3*nfreq,3))
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
        return np.zeros((1, 3*nfreq, 2*nin_max, 2*nin_max), dtype='complex')

problems = []
S = np.array([[read_single_S_list(filen, i ,j)[0] for j in range(nconf)] for i in range(nr)])
print 'WARNING: problems encountered in files', problems
#for prob in problems:
#    subprocess.call(["rm", data_direc+"Smat."+prob+".dat"])

def calc_q_mean(S_list, i, j):
    '''
    Calculates traces of Q operators of S matrices stored in S_list.
    If imaginary part of the traces are too high, their value is set
    to 0 and the number of configurations at this radius and frequency
    is decremented by 1.
    '''
    global dims
    global conf_counter
    global herm_counter
    
    q_mean_list=[]
    for w in range(nfreq):
        s = S_list[i,j][3*w:3*w+3]
        q_mean = np.trace(scat.calc_Q(s,dk)[0]) /(2*dims[w])
        q_mean_list.append(q_mean)            
    return q_mean_list

q_mean  = np.array([[calc_q_mean(S,i,j) for j in range(nconf)] for i in range(nr)])
means   = np.array([[np.mean(q_mean[i,:,j].real) for j in range(nfreq)] for i in range(nr)])
stddevs = np.array([[np.std(q_mean[i,:,j].real)  for j in range(nfreq)] for i in range(nr)])

def calc_qeig_mean(S_list, i, j):
    '''
    Calculates traces of Q operators of S matrices stored in S_list.
    If imaginary part of the traces are too high, their value is set
    to 0 and the number of configurations at this radius and frequency
    is decremented by 1.
    '''
    global dims
    global conf_counter
    global herm_counter
    
    q_mean_list=[]
    for w in range(nfreq):
        s = np.array(S_list[i,j,3*w:3*w+3])
        phases = np.zeros((3,), dtype='complex')
        #if (not np.all(s == 0.0)):
        for fnp in range(3):
            ###
            phases[fnp] = -imag_i * np.log(np.linalg.det(s[fnp]+10**(-32)))
        q_eig = (phases[2] - phases[0]) / (2*dk) / (2*dims[w])
        q_mean_list.append(q_eig)            
    return q_mean_list

q_mean  = np.array([[calc_qeig_mean(S,i,j) for j in range(nconf)] for i in range(nr)])
means   = np.array([[np.mean(q_mean[i,:,j].real) for j in range(nfreq)] for i in range(nr)])
stddevs = np.array([[np.std(q_mean[i,:,j].real)  for j in range(nfreq)] for i in range(nr)])

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
    plt.savefig(plot_direc+"AverageTime.Eigen."+filen+".%i.png" % indx, bbox_inches='tight')
    plt.clf()

        



  

