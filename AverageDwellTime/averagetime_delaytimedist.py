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
xi = -2./np.mean(ln_cond)*10.*L
v_eff = np.mean(np.sqrt(2.*energs[freq_range]))
N = int(np.mean(dims[freq_range]))
ls = 2.0/np.pi * xi / (N+1)
gamma = np.pi**2/4.0*ls
z = -0.5*np.mean(ln_cond)*10.*L
print "Gamma:", gamma
#print gamma*np.exp(z)
#print "Z", z
#print 1/z*np.log(10000./gamma)

ls1D = 0.5 * xi / (N+1)
gamma1D = 2.0 * ls
print "Gamma1D:", gamma1D


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
        mean_hist += hist[n] * time_pos * (time_range[-1]-time_range[0]) / bins
    mean_hist = mean_hist / (np.sum(hist) * (time_range[-1]-time_range[0]) / bins)
    return mean_hist

calc_bins_init = 8000
time_range = (max(np.min(q_eigv),10.**(-8)), np.max(q_eigv))
#time_range = (max(np.min(q_eigv),10.**(-8)), 175000.0)

calc_bins = calc_bins_init
mean_vals = np.mean(q_eigv)
calc_hist = np.histogram(q_eigv, bins=calc_bins, range=time_range, density=True)[0]
mean_hist = calc_mean_hist(time_range, calc_hist, calc_bins)
while(True):
    calc_bins_up = int(calc_bins*1.01)
    calc_hist_up = np.histogram(q_eigv, bins=calc_bins_up,
                                range=time_range, density=True)[0]
    mean_hist_up = calc_mean_hist(time_range, calc_hist_up, calc_bins_up)
    calc_bins_do = int(calc_bins*0.099)
    calc_hist_do = np.histogram(q_eigv, bins=calc_bins_do,
                                range=time_range, density=True)[0]
    mean_hist_do = calc_mean_hist(time_range, calc_hist_do, calc_bins_do)
    if ( abs(mean_vals - mean_hist) > abs(mean_vals - mean_hist_up) ):
        calc_bins = calc_bins_up
        calc_hist = calc_hist_up
        mean_hist = mean_hist_up
    elif ( abs(mean_vals - mean_hist) > abs(mean_vals - mean_hist_do) ):
        calc_bins = calc_bins_do
        calc_hist = calc_hist_do
        mean_hist = mean_hist_do
    else:
        break
    print "Bin Optimization:", calc_bins

calc_bins = 100

#fit_range = (10**(-8), 25000.0)
fit_range = (250.00, 9500.0)
#fit_range = time_range
fit_bins  = int((fit_range[-1] - fit_range[0]) / (time_range[-1] - time_range[0]) * calc_bins)
fit_hist = np.histogram(q_eigv, bins=fit_bins, range=fit_range, density=False)[0]
norm_hist = np.histogram(q_eigv, bins=calc_bins, range=time_range, density=False)[0]
fit_hist = fit_hist / (np.sum(norm_hist) * (time_range[-1]-time_range[0]) / calc_bins)
fitvals = np.linspace(fit_range[0], fit_range[-1], fit_bins)
print fit_bins

def fitfunc(tau, a):
    return 4. * a**2 / tau**3 * \
        np.exp(-2.*a/tau) * \
        (spspe.kn(0, 2*a/tau) + spspe.kn(1, 2*a/tau))

# print N*(N+1)
# def fitfunc(tau, N):
#     return 4. * (N*(N+1.)*gamma)**2 / tau**3 * \
#         np.exp(-2.*N*(N+1.)*gamma/tau) * \
#         (spspe.kn(0, 2*N*(N+1.)*gamma/tau) + spspe.kn(1, 2*N*(N+1.)*gamma/tau))

# def fitfunc(tau, a):
#     return a * 1./tau**2.

def fit_exp(tau, a, b, c):
    return a * c * np.exp(-b * tau**(1./1.))

# def fit_exp(tau, a, b, c):
#     return a *b* tau**c

fit_amp = spopt.curve_fit(fitfunc, fitvals, fit_hist)[0]
print "FIT PAR, GAMMA:", fit_amp, gamma

fit_a, fit_b, fit_c = spopt.curve_fit(fit_exp, fitvals, fit_hist)[0]
print "FIT PAR:", fit_a, fit_b, fit_c

mix_hist = np.zeros((calc_bins), dtype="float")

for n in range(calc_bins):
    time_pos = time_range[0] + (time_range[-1]-time_range[0]) / calc_bins * (n+0.5)

    if (time_pos <= fit_range[-1]):
        mix_hist[n] = calc_hist[n]
    else:
        mix_hist[n] = fitfunc(time_pos, fit_amp)
        mix_hist[n] = fit_exp(time_pos, fit_a, fit_b, fit_c)

obst_ratio = 0.94
Weyl = 10*L*10 * obst_ratio / (2*Pi)
print "Mean Vals:", np.mean(q_eigv)
print "Weyl:", Weyl
print "Mean Hist:", calc_mean_hist(time_range, calc_hist, calc_bins)
print "Mix  Hist:", calc_mean_hist(time_range, mix_hist, calc_bins)

dist_range = (time_range[0], 10000)
dist_bins  = int((dist_range[-1] - dist_range[0]) / (time_range[-1] - time_range[0]) * calc_bins)
dist_hist = np.zeros((dist_bins), dtype="float")
for n in range(dist_bins):
    time_pos = dist_range[0] + (dist_range[-1]-time_range[0]) / dist_bins * (n+0.5)
    dist_hist[n] = fitfunc(time_pos, 2.8*gamma)
distvals = np.linspace(dist_range[0], dist_range[-1], dist_bins)
print "Dist Hist:", calc_mean_hist(dist_range, dist_hist, dist_bins)

plt_range = time_range
plt_bins  = calc_bins
plt_hist  = np.histogram(q_eigv, bins=plt_bins, range=plt_range, density=True)[0]
pltvals = np.linspace(plt_range[0], plt_range[-1], plt_bins)

#plt.xlim(0.0,time_range[-1])
#plt.ylim(-1.0,0.2*np.mean(plt_hist))
plt.plot(pltvals, np.log10(plt_hist+10**(-10)), 'b-', label='orig_hist')
#plt.plot(fitvals, np.log10(fit_hist+10**(-10)), 'go', label='fit_hist')
#plt.plot(pltvals, np.log10(mix_hist+10**(-10)), 'g-', label='mix_hist')
#plt.plot(pltvals, np.log10(fitfunc(pltvals, fit_amp)), 'k-', label='fit')
#plt.plot(pltvals, np.log10(fitfunc(pltvals, 2.0*gamma)), 'r-', label='dist')
#plt.plot(distvals, np.log10(dist_hist+10**(-10)), 'y-', label='dist_hist')
plt.plot(pltvals, np.log10(mix_hist+10**(-10)), 'y-', label='mix_exp_hist')
plt.plot(pltvals, np.log10(fit_exp(pltvals,fit_a,fit_b,fit_c)+10**(-10)), 'g-', label='mix_exp')
plt.legend(bbox_to_anchor=(0.5, 1), loc=2, borderaxespad=0.)

plt.show()


  

