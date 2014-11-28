#!/usr/bin/env python

import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

from Scattering import scattering2D as scat


if (len(sys.argv)!=5):
    sys.exit("ABORT: parameters filen, # of radius variations, # of configurations,  # of frequency variations needed")

filen = str(sys.argv[1]) # namestring of calculation
nr    = int(sys.argv[2]) # number of radius variations
nconf = int(sys.argv[3]) # number of configurations to be averaged over
nfreq = int(sys.argv[4]) # number of frequency steps scanned over

data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/scatterdata/"

energs = scat.read_S(data_direc, filen+".0.0")[2]
kvals = map(np.sqrt, 2*energs)
dk = (kvals[2] - kvals[0]) / 2.0

nins = np.take(scat.read_S(data_direc, filen+".0.0")[1]/2, np.arange(1,3*nfreq,3))
nin_max = np.max(nins)

def read_single_S_list(filen, i, j):
    filen_i = filen+".%i.%i" % (i,j)
    return scat.read_S(data_direc, filen_i)

S = np.array([[read_single_S_list(filen, i ,j)[0] for j in range(nconf)] for i in range(nr)])

# t = []
# for i in range(nr):
#     for j in range(nfreq):
#         for k in range(nconf):
#             t.append(
#                 #np.take(np.take(S[i,k,1+3*j], range(dims[j]/2,dims[j]), axis=0), range(dims[j]/2), axis=1)
#                 S[i,k,1+3*j][dims[j]/2:dims[j],0:dims[j]/2]
#                 )
# # make numpy-array of t_amps, only valid for first three indices,
# # because then, the dimensions of the matrices change, generally
# t = np.array(t)
# t = np.reshape(t, (nr,nfreq,nconf))

t = np.zeros((nr, nfreq, nconf, nin_max, nin_max), dtype='complex')
for i in range(nr):
    for j in range(nfreq):
        for k in range(nconf):
            t[i,j,k,0:nins[j],0:nins[j]] = S[i,k,1+3*j][nins[j]:2*nins[j],0:nins[j]]

transmissions = np.array([[
            np.mean(np.absolute(t[i,j,:,:])) * nin_max**2 / nins[j]
            for j in range(nfreq)] for i in range(nr)])
#print transmissions

trans_mean = np.array([[
            np.mean(np.absolute(t[i,j,:,:])) * nin_max**2 / nins[j]**2
            for j in range(nfreq)] for i in range(nr)])
#print trans_mean

trans_msqr = np.array([[
            np.mean(np.absolute(t[i,j,:,:])**2) * nin_max**2 / nins[j]**2
            for j in range(nfreq)] for i in range(nr)])
#print trans_msqr

print (trans_msqr - trans_mean**2) / transmissions**2

#trans_msqr = np.vectorize(np.mean)(transmissions**2)

#print np.shape(transmissions)
#print np.shape(trans_mean)
#print transmissions[0,0,0]
#print trans_mean[0,0]

#print trans_mean
#print trans_msqr

#trans_vars = (trans_msqr/trans_mean)**2 - 1

#print trans_vars


def calc_q_mean(S_list, i, j):
    return [np.trace( scat.calc_Q(s,dk)[0] ) for s in np.reshape(S_list[i,j], (nfreq,3))]

q_mean  = np.array([[calc_q_mean(S,i,j) for j in range(nconf)] for i in range(nr)])
means   = np.array([[np.mean(q_mean[i,:,j]) for j in range(nfreq)] for i in range(nr)]).real
stddevs = np.array([[np.std(q_mean[i,:,j])  for j in range(nfreq)] for i in range(nr)]).real

plot_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/"

for mean, stddev, indx in zip(means, stddevs, range(nr)):

    plt.errorbar(range(nfreq), mean, stddev, fmt='D-r', label=r'$\tau_{tot}$')

    #Weyl = [10*2*10 * 0.96 / (2*Pi)] * nparam
    #plt.plot(freqs, Weyl, '-g', label=r"$\rm{Weyl-estimate}$", lw=2.0)

    plt.legend( bbox_to_anchor=(1.,1.) )
    plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
    plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    plt.ylabel(r'$\left\langle \rm{Tr}\left( Q \right) \right\rangle$')
    plt.savefig(plot_direc+"AverageTime."+filen+".%i.png" % indx, bbox_inches='tight')
    plt.clf()
 

#     plt.figure(2)

#     plt.errorbar(freqs, meankinds[1], stddevkinds[1], fmt='d-r', label=r'$\tau_{11}$')
#     plt.errorbar(freqs, meankinds[2], stddevkinds[2], fmt='d-g', label=r'$\tau_{t}$')
#     plt.errorbar(freqs, meankinds[3], stddevkinds[3], fmt='d-b', label=r'$\tau_{r}$')
#     polycoeffs = np.polyfit(freqs, meankinds[1], 5)
#     poly = np.poly1d(polycoeffs)
#     plt.plot(freqs, poly(freqs), '--y', label=r'$\tau_{11} \, \rm{fit}$', lw=2.0)
#     plt.legend( bbox_to_anchor=(1.,1.) )
#     plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
#     plt.ylabel(r'delay times')
#     plt.savefig("NLefttimes."+filen+".%i.png"%i, bbox_inches='tight')
#     plt.clf()

#     plt.figure(3)

#     plt.errorbar(freqs, meankinds[-1], stddevkinds[-1], fmt='d-r', label=r'$T$')
#     plt.errorbar(freqs, meankinds[-2], stddevkinds[-2], fmt='d-g', label=r'$R$')
#     plt.plot(freqs, np.array(meankinds[-1]) + np.array(meankinds[-2]), 'd--y', label=r'$T+R$')
#     #polycoeffs = np.polyfit(freqs, meankinds[1], 5)
#     #poly = np.poly1d(polycoeffs)
#     #plt.plot(freqs, poly(freqs), '--y', label=r'$\tau_{11} \, \rm{fit}$', lw=2.0)
#     plt.axis((freqs[0],freqs[-1], 0., 1.2))
#     plt.legend( bbox_to_anchor=(1.,1.) )
#     plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
#     plt.ylabel(r'$\rm{Transmittivity/Reflectivity}$')
#     plt.savefig("T+R."+filen+".%i.png"%i, bbox_inches='tight')
#     plt.clf()


# sp.call(["sh", "calcav.sh"])

    

        



  

