#!/usr/bin/env python

import numpy as np
import scipy as sp
import scipy.optimize as spopt

import matplotlib
matplotlib.use('Agg') # whatever, it makes it run via ssh

import matplotlib.pyplot as plt
import sys
from math import pi as Pi
import pickle

import warnings
import subprocess

from Utils import utils as ut
from Scattering import scattering2D as scat
from Calc import fitfuncs as ff
from Calc import transmission as trans
from Calc import localization as local
from Calc import phasetimes


if (len(sys.argv)!=6):
    sys.exit("ABORT: parameters filen, # of radius variations, # of configurations,  # of frequency variations, length of system [W] needed")

filen = str(sys.argv[1]) # namestring of calculation
nr    = int(sys.argv[2]) # number of radius variations
nconf = int(sys.argv[3]) # number of configurations to be averaged over
nfreq = int(sys.argv[4]) # number of frequency steps scanned over
L = float(sys.argv[5]) # length of system in units of width

data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/scatterdata/"

do_pickle = False

energs = scat.read_S(data_direc, filen+".0.0")[2]
if (len(energs)==0):
    kmin = 10.1 * np.pi / 10.
    kmax = 4.9 * np.pi / 10.
    dk = 0.5 * (kmin+kmax) * 10**(-7)
else:
    kvals = np.sqrt(2*energs)
    #dk = (kvals[-1] - kvals[-3]) / 2.0
    dk = (kvals[-1] - kvals[-2]) / 1.0
    #dk = 10.5*np.pi/10.*10.**(-7)

kvals = kvals[::3]
print "dk:", dk, 0.5*(10.1+14.9)*np.pi/10.*10.**(-10)
##### for av40
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

if (do_pickle):
    problems = []
    S = np.array([[read_single_S_list(filen, i ,j)[0] for j in range(nconf)] for i in range(nr)])
    print 'WARNING: problems encountered in files', problems

    t  = trans.calc_t (S, dims, (nr,nfreq,nconf,nin_max), fullout=True)
    r  = trans.calc_r (S, dims, (nr,nfreq,nconf,nin_max), fullout=True)
    tp = trans.calc_tp(S, dims, (nr,nfreq,nconf,nin_max), fullout=True)
    rp = trans.calc_rp(S, dims, (nr,nfreq,nconf,nin_max), fullout=True)

#     print "PICKLING scattering matrices"
#     pickle.dump( S, open(data_direc+"SMatrices."+filen+".p", "wb") )
# if (not do_pickle):
#     print "LOADING scattering matrices"
#     S = pickle.load( open(data_direc+"SMatrices."+filen+".p", "rb") )




conf_counter = nconf + np.zeros((nr, nfreq), dtype="int")
herm_counter = 0
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
                         
        if (abs(q_mean.imag) > 10**(-0)):
            q_mean = 0.0
            conf_counter[i,w] -= 1
            herm_counter += 1
            #print i,w,conf_counter[i,w]
        q_mean_list.append(q_mean)            
    return q_mean_list

if (do_pickle):
    q_mean  = np.array([[calc_q_mean(S,i,j) for j in range(nconf)] for i in range(nr)])

    qt_bar  = np.array([[[
                    np.trace(scat.calc_Q(t[i,3*w:3*w+3,j],dk)[0])/dims[w]
                    for j in range(nconf)] for w in range(nfreq)] for i in range(nr)])
    qr_bar  = np.array([[[
                    np.trace(scat.calc_Q(r[i,3*w:3*w+3,j],dk)[0])/dims[w]
                    for j in range(nconf)] for w in range(nfreq)] for i in range(nr)])
    qtp_bar = np.array([[[
                    np.trace(scat.calc_Q(tp[i,3*w:3*w+3,j],dk)[0])/dims[w]
                    for j in range(nconf)] for w in range(nfreq)] for i in range(nr)])
    qrp_bar = np.array([[[
                    np.trace(scat.calc_Q(rp[i,3*w:3*w+3,j],dk)[0])/dims[w]
                    for j in range(nconf)] for w in range(nfreq)] for i in range(nr)])

    print "PICKLING delay time data"
    pickle.dump( (q_mean, qt_bar, qr_bar, qtp_bar, qrp_bar), open(data_direc+"qMean."+filen+".p", "wb") )

    t  = t[:,np.arange(1,3*nfreq,3),:]
    tdt_eigvals = trans.calc_tdt_eigvals(t, (nr,nfreq,nconf,nin_max))
    print "PICKLING transmission data"
    pickle.dump( tdt_eigvals, open(data_direc+"Transmission."+filen+".p", "wb") )

if (not do_pickle):
    print "LOADING delay time data"
    q_mean, qt_bar, qr_bar, qtp_bar, qrp_bar = pickle.load( open(data_direc+"qMean."+filen+".p", "rb") )

    print "LOADING transmission data"
    tdt_eigvals = pickle.load( open(data_direc+"Transmission."+filen+".p", "rb") )


T_bar = np.array([[[
    np.sum(tdt_eigvals[i,j,k]) / dims[j]
    for k in range(nconf)] for j in range(nfreq)] for i in range(nr)])
R_bar = 1.0 - T_bar


qT_mean = np.mean( (qt_bar + qtp_bar) / (2.0*T_bar), axis=2).real
qR_mean = np.mean( (qr_bar + qrp_bar) / (2.0*R_bar), axis=2).real
qTR_mean = np.mean( 0.5 * (qt_bar + qtp_bar+ qr_bar + qrp_bar), axis=2).real


means   = np.array([[np.mean(q_mean[i,:,j].real) * nconf/max(conf_counter[i,j],1.0) for j in range(nfreq)] for i in range(nr)])
stddevs = np.array([[np.std(q_mean[i,:,j].real)  for j in range(nfreq)] for i in range(nr)])

local.check_single_channel(tdt_eigvals, 10, (nr,nfreq,nconf))

bins = np.array(50, dtype='float')
# range of which ln(G) values are considered for histogram
# must be < 0 otherwise error in localized fit
ln_cond_range = (-12,np.log(0.999))
# range out of which tdt-eigvals are considered for histogram
# must not be smaller than this otherwise error in localized fit
# for av25.1 best values for normalization in order to fit analytical
# formula are (0.15,0.85)
tdt_range = (0.15, 0.85)
tdt_full_range = (0.01, 0.99)

#mean_vars = local.calc_localization_meas_I(t, dims, (nr,nfreq,nin_max))
loc_var, cryst_hists = local.calc_localization_meas_II(tdt_eigvals, nr, bins)
ln_cond, lognorm_hists = local.calc_localization_meas_III(tdt_eigvals, (nr,nfreq,nconf), bins, ln_cond_range)
bins_full = bins * 19.0/(ln_cond_range[-1]-ln_cond_range[0])
ln_cond_full, lognorm_hists_full = local.calc_localization_meas_III(tdt_eigvals, (nr,nfreq,nconf), bins_full, (-16,3))

tdt_hists = local.calc_localization_meas_IV(tdt_eigvals, bins, tdt_range)
bins_full = bins * 1.0/(tdt_range[-1]-tdt_range[0])
#tdt_full_hists = local.calc_localization_meas_IV(tdt_eigvals, bins_full, tdt_full_range)
tdt_full_hists = local.calc_localization_meas_IV(tdt_eigvals, bins, tdt_full_range)
tdt_log_hists = local.calc_localization_meas_IV(np.log10(tdt_eigvals), bins, tdt_range, log=True)

    

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



### average time ###
for mean, stddev, indx in zip(means, stddevs, range(nr)):

    ### estimate derived from Weyl's Conjecture ###
    # number of open modes as integers
    dims_int = np.array(np.linspace(10.1, 14.9, nfreq, endpoint=True), dtype="int") 
    # number of open modes smoothened
    dims_float = np.array(np.linspace(10.1, 14.9, nfreq, endpoint=True), dtype="float")-0.5 

    # values for number of obstacles and their radii for av57 and av26
    Nobst = np.array([211, 13])
    radii = np.array([0.015*10, 0.060*10])

    A = 10.0 * L * 10.0 - Nobst[indx] * np.pi*radii**2
    B = 2.0 * 10.0 * (L - 1.0) + Nobst[indx] * 2.0 * np.pi * radii

    Weyl = 1.0 / (2.0 * dims_float) * (A[indx] * kvals - 0.5 * B[indx])
    Weyl_lead = 1.0 / (2.0 * (dims_float+0.5)) * (A[indx] * kvals)
    Weyl_lead_ana = np.pi * A[indx] / (2.0 * 10.0)

    print
    print "Total area:", 10.0 * L * 10.0
    print "Effective area:", A
    print "Number of obtacles:", Nobst
    print "Average Weyl time:", np.mean(Weyl)
    print "Difference of Weyl leading term and analytic exp.:", Weyl_lead[0] - Weyl_lead_ana
    print 

 

    plt.plot(range(nfreq), Weyl_lead, '--', color="#27c75f", label=r"$\rm{\tau_{Weyl-lead}}$ (smooth)", lw=3.0)
    plt.plot(range(nfreq), mean, '-k', label=r'$\tau_{tot}$', lw=3.0)
    plt.plot(range(nfreq), qT_mean[indx], '-b', label=r'$\tau_{T}$', lw=2.0)
    plt.plot(range(nfreq), qR_mean[indx], '-r', label=r'$\tau_{R}$', lw=2.0) 

    

    if (indx==0):
        #plt.semilogy()
        #plt.ylim(1e+1,1e+3)
        plt.xticks((0,24,48))
        plt.xlim(1., 48.)
        plt.ylim(20,70)
        plt.yticks((20,45))
        plt.plot(range(nfreq), Weyl, '--', color="#bf50bf", label=r"$\rm{\tau_{Weyl}}$ (smooth)", lw=3.0)
        plt.plot(range(nfreq), Weyl_lead, '--', color="#27c75f",
                 label=r"$\rm{\tau_{Weyl-lead}}$ (smooth)", lw=3.0)
        plt.plot(range(nfreq), mean, '-k', label=r'$\tau_{tot}$', lw=3.0)
        plt.plot(range(nfreq), qR_mean[indx], '-r', label=r'$\tau_{R}$', lw=2.0) 
        plt.savefig(plot_direc+"TimePartition."+filen+".%i_R.eps" % indx)
        plt.clf()

        plt.xticks((0,24,48))
        plt.xlim(1., 48.)
        plt.ylim(200,300)
        plt.yticks((200,300))
        plt.plot(range(nfreq), qT_mean[indx], '-b', label=r'$\tau_{T}$', lw=2.0)
        plt.savefig(plot_direc+"TimePartition."+filen+".%i_T.eps" % indx)
        plt.clf()
      
        
    if (indx==1):
        plt.xticks((0,24,48))
        plt.xlim(1., 48.)
        plt.ylim(20., 70.)
        plt.yticks((20,45,70))

        plt.savefig(plot_direc+"TimePartition."+filen+".%i.eps" % indx)
        plt.clf()

    plt.plot(range(nfreq), Weyl, '--', color="#27c75f", label=r"$\rm{\tau_{Weyl}}$ (smooth)", lw=3.0)
    plt.plot(range(nfreq), mean, '-k', label=r'$\tau_{tot}$', lw=3.0)
    plt.xticks((25,48))
    plt.yticks((32,34,36))
    plt.xlim(25., 48.)
    plt.ylim(32., 36.)
    plt.savefig(plot_direc+"TotalTime."+filen+".%i.eps" % indx)
    plt.clf()
 

### distribution for localized dynamics ###
for loc, hist, indx in zip(ln_cond, lognorm_hists, range(nr)):

    xvals = np.linspace(ln_cond_range[0],ln_cond_range[-1],bins).real
    bar_xvals =  np.linspace(ln_cond_range[0],ln_cond_range[-1],bins,endpoint=False).real
    dx = bar_xvals[1] - bar_xvals[0]
    mean_lnG = np.mean(ln_cond[indx]).real
    print "xi/L:", -2.0/mean_lnG
    def localized_log_amp_only_wrap(T, a): return a * ff.localized_log_amp_only(T, mean_lnG)

    ### renormalize histogram in considered range ###
    normalize = np.sum(hist)*dx
    if normalize==0: normalize = 1
    hist = 1.0 * hist / normalize


    try:
        #a = spopt.curve_fit(localized_log_amp_only_wrap, np.exp(xvals), hist.real)[0]
        #plt.plot(xvals, np.log10(localized_log_amp_only_wrap(np.exp(xvals), a)+1), label=r'localized dist')

        ### 16.2358: specific value for normalization of localized distribution with
        ### histogram-range (-12,0) and xi=0.43908357 calculated from av25.0.
        ### Value calculated with Mathematica file 'DistributionNormalizations.nb'.
        ### No fit involved here, comparison with analytical formula.
        plt.plot(xvals, -np.log10(localized_log_amp_only_wrap(np.exp(xvals), 1.0/16.2358)), color='blue', lw=4.0,
                 label=r'$-\rm{log}\left[\rho_{\rm{loc}}\right]$')
        pass
    except RuntimeError:
        print "WARNING: error in Loc fit III with radius #%i"%indx

    plt.bar(bar_xvals, -np.log10(hist.real+0.000001), width=dx, color='orange', label=r'$-\rm{log}\left[\rho\right]$') 

    #plt.legend( bbox_to_anchor=(-0.06,1.) )
    plt.xlim(-12, 0)  # range adjusted to av25.0
    plt.ylim(0.0,3.0) # range adjusted to av25.0
    plt.xticks((-12,-6,0))
    plt.yticks((0.0,1.5,3.0))
    #plt.title(r'$\rm{ln}\left[G_n\right]$-distribution for %i configurations at %i frequencies'%(nconf,nfreq))
    #plt.xlabel(r'$\rm{ln}\left[G_n\right]$')
    #plt.ylabel(r'$-\rm{log}\left[\rho(\rm{ln}\left[G_n\right])\right]$')
    plt.savefig(plot_direc+"Localization.III."+filen+".%i.eps" % indx, bbox_inches='tight')
    plt.clf()



### distribution for chaotic dynamics ###
for loc, hist, hist_full, hist_log, indx in zip(tdt_eigvals, tdt_hists, tdt_full_hists, tdt_log_hists, range(nr)):

     xvals = np.linspace(tdt_range[0],tdt_range[-1],bins)
     #xvals_full = np.linspace(0.0001,0.999,bins_full,endpoint=True)
     #bar_xvals =  np.linspace(0.0001,0.999,bins_full,endpoint=False)
     xvals_full = np.linspace(0.0001,0.999,bins,endpoint=True)
     bar_xvals =  np.linspace(0.0001,0.999,bins,endpoint=False)
     dx = bar_xvals[1] - bar_xvals[0]
     ### renormalize histogram in considered range ###
     normalize = np.sum(hist_full)*dx
     if normalize==0: normalize = 1
     hist = 1.0 * hist / normalize
     hist_full = 1.0 * hist_full / normalize


     try:
         #chaosamp = spopt.curve_fit(ff.chaotic, xvals, hist.real)[0]
         #plt.plot(xvals, ff.chaotic(xvals, chaosamp), label=r'chaotic dist')#(%.3f)'%cut_val)

         ### No fit involved here, Pi is analytical normalization constant.
         plt.plot(xvals_full, ff.chaotic(xvals_full, 1.0/np.pi), label=r'$\rho_{\rm{chaotic}}$', color='blue', lw=4.0)

         #plt.plot(xvals_full, ff.diffusive(xvals_full, 1.0/5.98645), label=r'$\rho_{\rm{chaotic}}$')
         pass
     except RuntimeError:
         print "WARNING: error in chaotic fit"



     plt.bar(bar_xvals, hist_full.real, width=dx, color='orange', label=r'$\rho$') 
     #plt.legend( bbox_to_anchor=(-0.06,1.) )
     #plt.title(r'$T_n$-distribution for %i configurations at %i frequencies'%(nconf,nfreq))
     #plt.xlabel(r'$T_n$')
     #plt.ylabel(r'$\rho(T_n)$')
     plt.xlim(0.0,1.0)
     plt.xticks((0.0,0.5,1.0))
     plt.ylim(0.0,3.0) # range adjusted to av25.1 
     plt.yticks((0.0,1.5,3.0))
     plt.savefig(plot_direc+"Localization.IV."+filen+".%i.eps" % indx, bbox_inches='tight')
     plt.clf()


    

        



  

