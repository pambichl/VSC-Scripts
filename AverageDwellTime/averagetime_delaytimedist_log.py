#!/usr/bin/env python

import numpy as np
import scipy as sp
import scipy.optimize as spopt
import scipy.special as spspe
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

import warnings
import subprocess

from Scattering import scattering2D as scat
from Utils import utils as ut

from Calc import transmission as trans

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

t = trans.calc_t(S, dims, (nr,nfreq,nconf,nin_max))

def calc_q_eigv(S_list, i, j, w):
    '''
    Calculates traces of Q operators of S matrices stored in S_list.
    If imaginary part of the traces are too high, their value is set
    to 0 and the number of configurations at this radius and frequency
    is decremented by 1.
    '''
    s = S_list[i,j][3*w:3*w+3]
    return ut.sort_eig(scat.calc_Q(s,dk)[0])[0].real

def calc_ln_cond(t_list, i, j, w):
    '''
    Calculates natural logarithm of total conductance
    for given transmission matrix.
    '''
    ln_cond = np.log(np.trace(np.dot(t_list[i,w,j].T.conj(), t_list[i,w,j])))
    return ln_cond.real


freq_range = range(9)
q_eigv = []
ln_cond = []
for i in range(nr):
    for j in range(nconf):
        for w in freq_range:
            q_eigv += list(calc_q_eigv(S,i,j,w))
            ln_cond.append(calc_ln_cond(t,i,j,w))

#np.save( "../../VSC2/AverageDwellTime-20130618/" + filen + '/ln_cond', ln_cond)
ln_cond = np.load("../../VSC2/AverageDwellTime-20130618/" + filen + '/ln_cond.npy')

q_eigv = np.array(q_eigv)
#np.save( "../../VSC2/AverageDwellTime-20130618/" + filen + '/delaytimedist', q_eigv)
q_eigv = np.load("../../VSC2/AverageDwellTime-20130618/" + filen + '/delaytimedist.npy')

def calc_mean_hist(time_range, hist, bins):
    '''
    Calculates mean value for given histogram.
    '''
    mean_hist = 0.0
    for n in range(bins):
        time_pos = time_range[0] + (time_range[-1]-time_range[0]) / bins * (n+0.5)
        mean_hist += hist[n] * 10**time_pos * (time_range[-1]-time_range[0]) / bins
    mean_hist = mean_hist / (np.sum(hist) * (time_range[-1]-time_range[0]) / bins)
    return mean_hist



print np.mean(q_eigv)
q_eigv = q_eigv[q_eigv > 0]
print np.mean(q_eigv)

# Zusammenspiel von bins, time_range, und fit_range!
# Ueberlegen, wann Erwartungswert ueberhaupt konvergiert!
# Zusatz durch Funktion nicht ueber Histogramm sondern direkt
# Integral ueber Funktion ausrechnen!

calc_bins = 2500
time_range = np.log10((np.min(q_eigv), 10000000000*np.max(q_eigv)))
timevals = np.linspace(time_range[0], time_range[-1], calc_bins)

calc_hist = np.histogram(np.log10(q_eigv), bins=calc_bins, range=time_range, density=True)[0]
mean_hist = calc_mean_hist(time_range, calc_hist, calc_bins)
print mean_hist

def fit_func(tau, a, b, c, d, e):
    #return a * (tau)**e * np.exp(-b * (tau)**d)
    return a * np.exp(-b * (tau)**d)

def Schomerus_func(tau, a):
    return tau * np.log(10) * 4. * a**2 / tau**3 * \
        np.exp(-2.*a/tau) * \
        (spspe.kn(0, 2*a/tau) + spspe.kn(1, 2*a/tau))


fit_range = (2.0, 4.0)
fit_bins  = int((fit_range[-1] - fit_range[0]) / (timevals[-1] - timevals[0]) * calc_bins)
fit_hist = np.histogram(np.log10(q_eigv), bins=fit_bins, range=fit_range, density=False)[0]
norm_hist = np.histogram(np.log10(q_eigv), bins=calc_bins, range=time_range, density=False)[0]
fit_hist = fit_hist / (np.sum(norm_hist) * (time_range[-1]-time_range[0]) / calc_bins)
fitvals = np.linspace(fit_range[0], fit_range[-1], fit_bins)

a, b, c, d, e = spopt.curve_fit(fit_func, fitvals, fit_hist, maxfev = 100000)[0]
print a, b, c, d, e

mix_hist = np.zeros((calc_bins), dtype="float")
for n in range(calc_bins):
    time_pos = time_range[0] + (time_range[-1]-time_range[0]) / calc_bins * (n+0.5)
    if (time_pos <= fit_range[-1]):
        mix_hist[n] = calc_hist[n]
    else:
        mix_hist[n] = fit_func(time_pos, a, b, c, d, e)

mean_hist = calc_mean_hist(time_range, mix_hist, calc_bins)
print mean_hist

plt.ylim(0.0,2.0)
plt.plot(timevals, calc_hist, 'b-', label='orig_hist')
plt.plot(fitvals, fit_hist, 'y-', label='fit_hist')
plt.plot(timevals, fit_func(timevals, a, b, c, d, e), 'g-', label='fit')
plt.plot(timevals, Schomerus_func(timevals, 0.68), 'k--', label='Schomerus')
plt.plot(timevals, mix_hist, 'r-', label='mix_hist')
plt.legend(bbox_to_anchor=(0.5, 1), loc=2, borderaxespad=0.)

plt.show()


  

