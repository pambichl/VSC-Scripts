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

if (len(sys.argv)!=5):
    sys.exit("ABORT: parameters filen, # of radius to consider, # of configurations, # of frequency steps needed.")

filen = str(sys.argv[1]) # namestring of calculation
nr    = int(sys.argv[2]) # radius variation to consider
nconf = int(sys.argv[3]) # number of configurations to be averaged over
nfreq = int(sys.argv[4]) # number of frequencies scanned

data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"
plot_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"


lead_width = 10.
kmin = 10.5 * np.pi / lead_width
kmax = 10.5 * np.pi / lead_width
dk = 0.5 * (kmax+kmin) * 10.**(-10)
k_range = np.linspace(kmin, kmax, nfreq)
E_range = 0.5 * k_range**2
dims = np.array(np.floor(k_range * lead_width / np.pi), dtype='float') # float or int??
#dims = np.array(k_range * lead_width / np.pi, dtype='float')

def read_single_EVal_list(filen, i, j, w):
    '''
    Read in list of complex eigenvalues of
    effective Hamiltonian from file.
    '''
    global problems
    
    filen_i = "EVal_C++." + filen + ".%i.%i.%04i.dat" % (i,j,w)
    location = data_direc + "resonancedata/" + filen_i
    print "reading file: " + filen_i
    #try:
    data = np.loadtxt(location, dtype='float', unpack=True)
    RealSort = np.argsort(data[1])
    data = data[:,RealSort]
    NVals = int(1.0*np.shape(data)[1])
    #return data
    print data[:,0]
    return data[:,:NVals]

    # except:
    #     if 'problems' in globals():
    #         problems.append(filen_i) 
    #     print "WARNING: problem reading "+filen_i
    #     return 0.0

problems = []
data = np.array([[read_single_EVal_list(filen, nr, j, w) for j in range(nconf)] for w in range(nfreq)])
print 'WARNING: problems encountered in files', problems

Reals = data[:,:,1,:]
Imags = data[:,:,2,:]
print np.shape(Reals)

pole_mask = (Imags >= 0.0)
Reals[pole_mask] = 1000.*w
Imags[pole_mask] = -10.**(-32)
EVals = Reals + imag_i*Imags

Tpref = 2.0
Ipref = 0.5 * 1./np.pi
Tpref = 2.0
Ipref = 1.0
Imags = Ipref * Imags

#print 1./(E_range[-1]-np.max(Reals))**2, 1./(E_range[0]-np.min(Reals))**2

tau_W_list = np.array(
    [k_range[w]/dims[w] * np.sum( -Imags[w] / ((E_range[w] - Reals[w])**2 + (Imags[w])**2), axis=1) for w in range(nfreq)])
# tau_W_list = np.array(
#     [Tpref * k_range[w] * 1./(2*dims[w]) * np.sum( -Imags[w] / ((E_range[w] - Reals[w])**2 + Imags[w]**2), axis=1) for w in range(nfreq)])
tau_W = np.mean(tau_W_list, axis=1)
stddev = np.std(tau_W_list, axis=1)
#print "no integration:", tau_W

Ediff = np.max(Reals) - np.min(Reals)
rel = 0.075
E_intvl = np.array([E_range*(1.0-rel), E_range*(1.0+rel)]).T
N_pos = 1000
dE = (E_intvl[:,1]-E_intvl[:,0]) / N_pos

# print "dE:", dE
# print "k-range:", np.array([np.sqrt(2.*E_intvl[:,0]), k_range, np.sqrt(2.*E_intvl[:,1])]).T

# mean_tau_W_list = np.zeros((nfreq, nconf), dtype='float')
# mean_tau_W_list_pref = np.zeros((nfreq, nconf), dtype='float')
# for w in range(nfreq):
#     for p in range(N_pos):
#         pos = E_intvl[w,0] + p*dE[w]
#         # numerical integration
#         mean_tau_W_list += 1./N_pos * \
#         np.sqrt(2.*pos)/dims[w] * \
#         np.sum( -Imags[w] / ((pos - Reals[w])**2 + Imags[w]**2), axis=1)
#         # numerical integration with approximate prefactor
#         mean_tau_W_list_pref += 1./N_pos * \
#         k_range[w]/dims[w] * \
#         np.sum( -Imags[w] / ((pos - Reals[w])**2 + Imags[w]**2), axis=1)
# mean_tau_W = np.mean(mean_tau_W_list, axis=1)
# print "numerical:", mean_tau_W
# mean_tau_W_pref = np.mean(mean_tau_W_list_pref, axis=1)
# print "numerical (approx. pref.):", mean_tau_W_pref

# analytical integration (error since k=k(E), but probably small
mean_tau_W_list = np.zeros((nfreq, nconf), dtype='float')
for w in range(nfreq):
    mean_tau_W_list[w] = 2.0*k_range[w]/(2*dims[w]) * 1./(E_intvl[w,1]-E_intvl[w,0]) * np.sum(
    np.arctan((E_intvl[w,1]-Reals[w])/(-0.5/(2*np.pi)*Imags[w])) - np.arctan((E_intvl[w,0]-Reals[w])/(-0.5/(2*np.pi)*Imags[w])),
    axis=1)
mean_tau_W = np.mean(mean_tau_W_list, axis=1)
Int_stddev = np.std(mean_tau_W_list, axis=1)
#print "analytical:", mean_tau_W

# def calculate_freq_mean(w, w_intvl):
#     global k_range, Reals, Imags
#     mean_tau_W_list = k_range[w]/dims[w] * 1.0/(w_intvl[1]-w_intvl[0]) * np.sum(
#         np.arctan((w_intvl[1]-Reals[w])/-Imags[w]) - np.arctan((w_intvl[0]-Reals[w])/-Imags[w]),
#         axis=1)
#     mean_tau_W = np.mean(mean_tau_W_list, axis=0)
#     return mean_tau_W

def calculate_freq_mean(w, w_intvl):
    global k_range,  dims, Reals, Imags
    mean_tau_W_list = k_range[w]/dims[w] * 1.0/(w_intvl[1]-w_intvl[0]) * np.sum(
        np.arctan((w_intvl[1]-Reals[w])/(-Imags[w])) - np.arctan((w_intvl[0]-Reals[w])/(-Imags[w])),
        axis=1)
    mean_tau_W = np.mean(mean_tau_W_list, axis=0)
    return mean_tau_W

rel_range = np.linspace(0.001, 0.20, 100)
for w in range(nfreq):
    E_pos = E_range[w]
    #print w, E_pos
    w_range = E_pos * np.array([1.0-rel_range, 1.0+rel_range]).T # wieso hier transponiert?
    Int_vec = np.array([calculate_freq_mean(w, w_range[i]) for i in range(len(rel_range))])
    plt.plot(rel_range, Int_vec, '-x', label="Eigenvalues" )
plt.savefig(plot_direc+"Res_Int."+filen+".png", bbox_inches='tight')
plt.clf()

#print "slgjsalgjslg", dims


try:
    obstacles = np.loadtxt(plot_direc + 'potential.dat', usecols=(1,))
    obst_mask = abs(obstacles) == 0
    obst_ratio = float(np.shape(obstacles[obst_mask])[0]) / float(np.shape(obstacles)[0])
    print "filling factor of obstacles: %5f" % (1.0-obst_ratio)
except:
    obst_ratio = 0.94
    print "WARNING: automatically assuming filling fraction to be 0.06"


plt.errorbar(range(nfreq), mean_tau_W, Int_stddev, fmt='D-b', label=r'$\tau_{tot}$')
plt.errorbar(range(nfreq), tau_W, stddev, fmt='D-r', label=r'$\tau_{tot}$')

L = 2.5

# Weyl estimate for dwell time: A / 2Pi = W * L * (1-fill) / 2Pi
Weyl = [10*L*10 * obst_ratio / (2*Pi)]*nfreq
plt.plot(range(nfreq), Weyl, '--g', label=r"$\rm{Weyl-estimate}$", lw=2.0)

mean_mean = np.mean(mean_tau_W)
plt.plot(range(nfreq), [mean_mean]*nfreq, '--', color='lightblue', label=r"$\rm{mean}$", lw=2.0)
print "Weyl estimate, Int-mean-value:", Weyl[0], mean_mean
mean_mean = np.mean(tau_W)
plt.plot(range(nfreq), [mean_mean]*nfreq, '--', color='orange', label=r"$\rm{mean}$", lw=2.0)
print "Weyl estimate, mean-value:", Weyl[0], mean_mean
print tau_W

plt.legend( bbox_to_anchor=(1.,1.) )
plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
plt.ylabel(r'$\left\langle \tau_W \right\rangle$')
plt.ylim(0., 2*Weyl[0])
plt.savefig(plot_direc+"AverageTime.Res."+filen+".%i.png" % nr, bbox_inches='tight')
plt.clf()


Hits = np.logical_and(Reals + 1*Imags <= w, Reals - 1*Imags >= w) # Imags < 0

indx = np.argsort(Imags.flatten())
EVals = EVals.flatten()[indx]
Hits = Hits.flatten()[indx]
HitEVals = EVals[Hits]

# print "EVals"
# print EVals[:10]
# print EVals[-10:]
# print "HitEVals"
# print HitEVals[:10]
# print HitEVals[-10:]
# print "k, dk", k, dk
# print "w, dw", w, dw

# vals_range = -np.log10(-np.array([np.min(EVals.imag), np.max(EVals.imag)]))
# hitvals_range = -np.log10(-np.array([np.min(HitEVals.imag), np.max(HitEVals.imag)]))
                                   

# bins = 100

# hist = np.histogram(-np.log10(-EVals.imag), bins=bins, range=vals_range, density=True)[0]
# hithist = np.histogram(-np.log10(-HitEVals.imag), bins=bins, range=vals_range, density=True)[0]

# ### Plots ###

# xvals = np.linspace(vals_range[0], vals_range[1], bins)

# plt.plot(xvals, hist, '-g', label=r"Imags", lw=1.5)
# plt.plot(xvals, hithist, '-y', label=r"HitImags", lw=1.5)
# #plt.plot(xvals, hithist/hist, '-r', label=r"Hit", lw=1.5)

# #plt.xlim(0.0,1.0)

# bins = 100
# x_pos = np.linspace(np.min(Reals), np.max(Reals), bins)
# hist = np.histogram(Reals.flatten(), bins=bins, density=True)[0]
# #plt.plot(x_pos, hist)

for w in range(nfreq):
    plt.plot( Reals[w].flatten(), np.log10(-Imags[w].flatten()), 'x', label="Eigenvalues" )
    plt.vlines(E_range[w], -30, 5, 'y')
    #plt.xlim((E_intvl[0,0], E_intvl[-1,-1]))

#plt.legend( bbox_to_anchor=(1.,1.) )
plt.title(r'bla')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig(plot_direc+"Resonances."+filen+".png", bbox_inches='tight')
plt.clf()

print "ATTENTION to dim (freq-range), prefactors, etc.?"
