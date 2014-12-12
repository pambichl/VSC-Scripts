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
from Utils import utils as ut

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
plot_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"

try:
    obstacles = np.loadtxt(plot_direc + 'potential.dat', usecols=(1,))
    obst_mask = abs(obstacles) == 0
    obst_ratio = float(np.shape(obstacles[obst_mask])[0]) / float(np.shape(obstacles)[0])
    print "filling factor of obstacles: %5f" % (1.0-obst_ratio)
except:
    obst_ratio = 0.94
    print "WARNING: automatically assuming filling fraction to be 0.06"

# Weyl estimate for dwell time: A / 2Pi = W * L * (1-fill) / 2Pi
Weyl = [10*L*10 * obst_ratio / (2*Pi)]*nfreq

energs = scat.read_S(data_direc, filen+".0.0")[2]
if (len(energs)==0):
    kmin = 3.1 * np.pi / 10.
    kmax = 4.9 * np.pi / 10.
    dk = 0.5 * (kmin+kmax) * 10**(-7)
else:
    kvals = np.sqrt(2*energs)
    dk = (kvals[-1] - kvals[-3]) / 2.0

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

r_list = trans.calc_r(S, dims, (nr,nfreq,nconf,nin_max), fullout=True)
t_list = trans.calc_t(S, dims, (nr,nfreq,nconf,nin_max), fullout=True)

def calc_q_mean(S_list, i, j):
    '''
    '''
    global dims
    global Weyl
    global r_list, t_list
    
    q_mean_list=[]
    for w in range(nfreq):
        s = S_list[i,j][3*w:3*w+3]
        Q = scat.calc_Q(s,dk)[0]
        qlist = ut.sort_eig(Q)[0]

        r = r_list[i,3*w:3*w+3,j,0:dims[w],0:dims[w]]
        t = t_list[i,3*w:3*w+3,j,0:dims[w],0:dims[w]]
        
        Q_eigstates = (ut.sort_eig(Q)[1].T)
        
        # qmax = np.sum(
        #     np.absolute(np.dot(r[1], state_max))**2 * \
        #         ( np.angle(np.dot(r[2], state_max)) - np.angle(np.dot(r[0], state_max)) ) / (2*dk) +\
        #         np.absolute(np.dot(t[1], state_max))**2 * \
        #         ( np.angle(np.dot(t[2], state_max)) - np.angle(np.dot(t[0], state_max)) ) / (2*dk)
        #     )

        for n, state in enumerate(Q_eigstates):
            qman = np.sum(
                np.absolute(np.dot(s[1], state))**2 * \
                    ( np.angle(np.dot(s[2], state)) - np.angle(np.dot(s[0], state)) ) / (2*dk)
                )
            #print np.shape(t[2]), np.shape(t[0]), np.shape(state), np.shape(state[0:dims[w]]), 
            tran = ( np.linalg.norm(np.dot(t[0], state[0:dims[w]]))**2,
                     np.linalg.norm(np.dot(t[1], state[0:dims[w]]))**2,
                     np.linalg.norm(np.dot(t[2], state[0:dims[w]]))**2 )
            refl = ( np.linalg.norm(np.dot(r[2], state[0:dims[w]]))**2, np.linalg.norm(np.dot(r[0], state[0:dims[w]]))**2 )
            if ( abs((qman-qlist[n])/qlist[n]) > 0.01 ):
                print "TIME", i, j, w, n, qman, qlist[n]
                qlist[n] = abs(qman)
            if ( abs((tran[1]-tran[0])/tran[0]) > 0.001 or
                 abs((tran[2]-tran[1])/tran[1]) > 0.001):
                print "TRANS", i, j, w, n, tran[1], tran[0], qman

        qmax = 0.0*qlist[-1]
        #qmax = 100000.0*Weyl[0]
        sdiff = np.array([s[0] * np.exp(imag_i*qmax*dk), s[1], s[2] * np.exp(-imag_i*qmax*dk)])
        Qdiff = scat.calc_Q(sdiff,dk)[0]
        #qlist[-1] = ut.sort_eig(Qdiff)[0][-1] + qmax
        #qlist = ut.sort_eig(Qdiff)[0] + qmax
        q_mean = np.sum(qlist)
        q_mean_list.append(q_mean)            
    return q_mean_list

q_mean  = np.array([[calc_q_mean(S,i,j) for j in range(nconf)] for i in range(nr)])
means   = np.array([[np.mean(q_mean[i,:,j]) / (2*dims[j]) for j in range(nfreq)] for i in range(nr)]).real
stddevs = np.array([[np.std(q_mean[i,:,j].real)  for j in range(nfreq)] for i in range(nr)])


tdt_eigvals = trans.calc_tdt_eigvals(t_list, (nr,nfreq,nconf,nin_max))


### Plots ###

for mean, stddev, indx in zip(means, stddevs, range(nr)):

    plt.errorbar(range(nfreq), mean, stddev, fmt='D-r', label=r'$\tau_{tot}$')

    plt.plot(range(nfreq), Weyl, '--g', label=r"$\rm{Weyl-estimate}$", lw=2.0)

    mean_mean = np.mean(mean)
    plt.plot(range(nfreq), [mean_mean]*nfreq, '--', color='orange', label=r"$\rm{mean}$", lw=2.0)
    print "Weyl estimate, mean-value:", Weyl[0], mean_mean

    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'$\left\langle \rm{Tr}\left( Q \right) \right\rangle$')
    plt.ylim(0., 2*Weyl[0])
    plt.savefig(plot_direc+"AverageTime_Diff."+filen+".%i.png" % indx, bbox_inches='tight')
    plt.clf()
 

        



  

