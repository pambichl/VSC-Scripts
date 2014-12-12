# S0 = scat.read_S(data_direc, filen+"_clean")[0]
# q0 = calc_q_mean( np.array([[S0]]), 0, 0)
# print np.shape(q0)


def theta_function(x):
    '''Heaviside-Theta function.'''
    if x == 0:
        return 0.5
    elif x < 0:
        return 0.0
    else:
        return 1.0

def ballistic_time(k,n):
    '''Analytical expression for n-th delay time
    for perfect (empty) waveguide.'''
    global L
    ky = n*np.pi/10.0
    time = L * 10.0 * k / np.sqrt(np.abs(k**2 - ky**2)) * theta_function(k-ky)
    return time

def ballistic_TrQ(k):
    '''Analytical expression for trace of Q-operator
    for perfect (empty) waveguide.'''
    Tr = 0.0
    for n in np.arange(1,15):
        Tr += ballistic_time(k,n)
    return 2.0*Tr


    ### values of actual calculations W = 10, L = 2.5 *W, nx = L/dx, ny = 14.9*10, pphwl = 10
    nx = 375; ny = 149; dx = 0.0666666666; ATot = nx*(ny+1)*dx**2
    AObst = ATot * (1.0-obst_ratio)
    ARem  = ATot * obst_ratio
    radii = np.array((0.015*10., 0.060*10.), dtype="float") # used sizes of obstacles r=0.015*W, r=0.060*W
    NObst = np.array(0.06*ATot / (np.pi*radii**2.), dtype="int") # number of obstacles, obst_ratio = 0.06 analytically
    BObst = NObst * 2. * radii * np.pi # surface of obstacles
    ### estimate for prefactor for correction of surface of obstacles
    #BPref = np.array( ( (12*0.0666666 / (2*radii[0]*np.pi)), 1.0 ) )
    BEff = BObst[indx] + 2.0 * nx * dx * L - 2.0 * (ny+1) * dx
    
    Weyl = (ARem * kvals / (2.0)
            - 1.0 * BObst[indx]/(2.0)
            - 2.*10.*L/(2.0)
            + 2.*10./(2.0))

    DOS = ARem * kvals / (2.0 * np.pi) - BEff / (4.0 * np.pi)
    Time = 2.0 * np.pi * DOS


    
    Weyl_step = Time/(2.0*dims_int)
    Weyl_cont = Time/(2.0*dims_float)
    #Weyl_lead = ARem *  np.pi / (2.0*10) * np.ones((nfreq))
    Weyl_lead = ARem *  kvals / (2.0*dims_float)

    Blanco = 4.0 * ARem / (2.0*10)
    Weyl_cont = Weyl


    #plt.plot(range(nfreq), Weyl_step, '--b', label=r"$\rm{\tau_{Weyl}}$ (step)", color='blue', lw=2.5)
    #plt.plot(range(nfreq), Weyl_lead, '--', label=r"$\rm{\tau_{Weyl}}$ (1st term)", color='orange', lw=2.5)
    ###   ###



    mean_mean = np.mean(mean)

    #plt.plot(range(nfreq), [mean_mean]*nfreq, '--', color='orange', label=r"$\rm{mean}$", lw=2.0)
    #print "Weyl estimate, mean-value:", Weyl[0], mean_mean


    #plt.plot(range(nfreq), np.ones((nfreq,))*mean_mean, '-', label="Mean", lw=1.5)



    # plt.errorbar(range(nfreq), mean, stddev, fmt='D-g', label=r'$\tau_{tot}$')
    # plt.plot(range(nfreq), Weyl_cont, '--g', label=r"$\rm{\tau_{Weyl}}$ (smooth)", color='red', lw=2.5)

    # plt.legend( bbox_to_anchor=(1.,1.) )
    # plt.title(r'average delay time for %i configurations at each $\omega$' % nconf)
    # plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    # plt.ylabel(r'$\left\langle \rm{Tr}\left( Q \right) \right /M \rangle$')
    # plt.xlim(0., 50.) # range adjusted to av26 and av57
    # plt.ylim(0., 70.) # range adjusted to av26 and av57
    # plt.savefig(plot_direc+"AverageTime."+filen+".%i.png" % indx, bbox_inches='tight')
    # plt.savefig(plot_direc+"AverageTime."+filen+".%i.eps" % indx, bbox_inches='tight')
    # plt.clf()





    #plt.legend( bbox_to_anchor=(1.,1.) )
    #plt.title(r'delay time partitioning for %i configurations at each $\omega$' % nconf)
    #plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    #plt.savefig(plot_direc+"TimePartition."+filen+".%i.png" % indx, bbox_inches='tight')
    #plt.savefig(plot_direc+"TimePartition."+filen+".%i.eps" % indx, bbox_inches='tight')



# for mean, stddev, indx in zip(meansIm, stddevsIm, range(nr)):

#     plt.errorbar(range(nfreq), mean, stddev, fmt='D-r', label=r'$\rm{Im}(\tau_{tot})$')

#     mean_mean = np.mean(mean)
#     plt.plot(range(nfreq), [mean_mean]*nfreq, '--', color='orange', label=r"$\rm{mean}$", lw=2.0)

#     plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#     #plt.ylim(-10**(-5),10**(-5))
#     plt.legend( bbox_to_anchor=(1.,1.) )
#     plt.title(r'Im of average delay time for %i conf. per $\omega$' % nconf)
#     plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
#     plt.ylabel(r'$\rm{Im}(\left\langle \rm{Tr}\left( Q \right) \right\rangle$)')
#     plt.savefig(plot_direc+"AverageTimeIm."+filen+".%i.png" % indx, bbox_inches='tight')
#     plt.clf()

# for loc, indx in zip(mean_vars, range(nr)):

#     plt.plot(range(nfreq), loc, 'D-r', label=r'$\left\langle s_{ab} \right\rangle$')
#     plt.plot(range(nfreq), [7./3.]*nfreq, '--g', label=r'$7/3$')

#     plt.legend( bbox_to_anchor=(1.,1.) )
#     plt.title(r'$\left\langle s_{ab} \right\rangle$ measure for localization for %i configurations at each $\omega$' % nconf)
#     plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
#     plt.ylabel(r'$\left\langle s_{ab} \right\rangle$')
#     plt.savefig(plot_direc+"Localization.I."+filen+".%i.eps" % indx, bbox_inches='tight')
#     plt.clf()


# for loc, hist, indx in zip(loc_var, cryst_hists, range(nr)):

#     plt.plot(np.linspace(np.min(loc.real),np.max(loc.real),bins), hist.real, label=r'$N(x_n)$')

#     plt.legend( bbox_to_anchor=(1.,1.) )
#     plt.title(r'$x_n$-distribution measure for localization for %i configurations at %i frequencies'%(nconf,nfreq))
#     plt.xlabel(r'$x_{n}$')
#     plt.ylabel(r'$N(x_n)$')
#     plt.savefig(plot_direc+"Localization.II."+filen+".%i.png" % indx, bbox_inches='tight')
#     plt.clf()



    # try:
    #     x0, s, a = spopt.curve_fit(ff.gauss, xvals, hist.real)[0]
    #     plt.plot(xvals, np.log10(ff.gauss(xvals, x0, s, a)), label=r'Gauss-fit')
    #     print "xi Gauss III [L]: %5f" % (-2./x0)
    # except RuntimeError:
    #     print "WARNING: error in Gauss fit III with radius #%i"%indx
    # try:
    #     a, xi = spopt.curve_fit(ff.localized_log, np.exp(xvals), hist.real)[0]
    #     plt.plot(xvals, np.log10(ff.localized_log(np.exp(xvals), a, xi)), label=r'Loc-fit')
    #     print "fit xi III [L]: %5f" % xi
    # except RuntimeError:
    #     print "WARNING: error in Loc fit III with radius #%i"%indx



#     loc = np.sort(loc.flatten())
#     tot_len = np.shape(loc)[0]

#     cut_val = 0.0

#     zero_mask = loc < cut_val
#     cut_len_zero = np.shape(loc[zero_mask])[0] # number of values below threshold
#     N_vals = 0
#     bin_cutoff_zero = 0
#     while(N_vals < cut_len_zero):
#         N_vals += hist[bin_cutoff_zero]
#         bin_cutoff_zero += 1
#     zero_mask = loc > (1-cut_val)
#     cut_len_one = np.shape(loc[zero_mask].flatten())[0] # number of values below threshold
#     N_vals = 0
#     bin_cutoff_one = 0
#     while(N_vals < cut_len_one):
#         N_vals += hist[-(bin_cutoff_one+1)]
#         bin_cutoff_one += 1

#     print "number of bins cut away at T=0 and T=1:", bin_cutoff_zero, bin_cutoff_one



#     x_logvals = np.linspace(np.log10(xvals[0]), np.log10(xvals[-1]), bins)
#     fit_inds = np.arange(bin_cutoff_zero, bins-(bin_cutoff_one), dtype='int')
#     mean_lnG = np.mean(ln_cond[indx]).real
#     def localized_amp_only_wrap(T, a): return a * ff.localized_amp_only(T, mean_lnG)
#     mean_lnT = np.mean(np.log(loc[cut_len_zero:tot_len-cut_len_one])).real
#     #mean_lnT = np.mean(-np.log(loc)).real

#     hist = np.array(hist, dtype='float')
#     hist_log = np.array(hist_log, dtype='float')

#     zero_hist_mask = hist_log == 0.0
#     hist_log[zero_hist_mask] = 0.001
#     zero_hist_mask = hist == 0.0
#     hist[zero_hist_mask] = 0.001

#     try:
#         diffamp = spopt.curve_fit(ff.diffusive, xvals[fit_inds], hist.real[fit_inds])[0]
#         plt.plot(xvals, ff.diffusive(xvals, diffamp), label=r'diffusive-fit(%.3f)'%cut_val)
#     except RuntimeError:
#         print "WARNING: error in diffusive fit"
     # try:
     #     chaosamp = spopt.curve_fit(ff.chaotic, xvals[fit_inds], hist.real[fit_inds])[0]
     #     plt.plot(xvals, ff.chaotic(xvals, chaosamp), label=r'chaotic-fit(%.3f)'%cut_val)
     # except RuntimeError:
     #     print "WARNING: error in chaotic fit"
#     try:
#         locamp, xi = spopt.curve_fit(ff.localized, xvals[fit_inds], hist.real[fit_inds])[0]
#         plt.plot(xvals, ff.localized(xvals, locamp, xi), label=r'single-channel-fit(%.3f)'%cut_val)
#         print "fit xi IV, eff. xi (G) [L]: %5f, %5f" % (xi, -2./mean_lnG)
#     except RuntimeError:
#         print "WARNING: error in localized fit"
#     try:
#         locamp_only = spopt.curve_fit(localized_amp_only_wrap, xvals[fit_inds], hist.real[fit_inds])[0]
#         plt.plot(xvals, ff.localized_amp_only(xvals, locamp_only), label=r'single-channel-fit(amp. only, %.3f)'%cut_val)
#     except:
#         print "WARNING: error in localized (amplitude only) fit"

     #plt.plot(xvals, hist.real, label=r'$N(T_n)$')
     #plt.plot(xvals_full, hist_full.real, label=r'$N(T_n)$')#



#     try:
#         loc_log_amp_only = spopt.curve_fit(localized_amp_only_wrap, 10**x_logvals, hist_log.real)[0]
#         loc_log_amp, xi_log = spopt.curve_fit(ff.localized_log, 10**x_logvals, hist_log.real)[0]
#         x0, s, a = spopt.curve_fit(ff.gauss, x_logvals, hist_log.real)[0]
#         plt.plot(x_logvals, np.log10(localized_amp_only_wrap(10.**x_logvals, loc_log_amp_only)), label=r'loc-fit(amp-only)')
#         plt.plot(x_logvals, np.log10(ff.localized_log(10.**x_logvals, loc_log_amp, xi_log)), label=r'loc-fit')
#         plt.plot(x_logvals, np.log10(ff.gauss(x_logvals, x0, s, a)), label=r'gauss-fit')
#         print "fit xi_log IV, xi_log Gauss [L]: %5f, %5f" % (xi_log, -2./x0)
#     except RuntimeError:
#         print "WARNING: error in localized (log, amplitude only) fit radius #%i"%indx

#     plt.plot(x_logvals, np.log10(hist_log.real), label=r'$\rm{log}(N(\rm{log}(T_n)))$')
#     plt.legend( bbox_to_anchor=(1.3,1.) )
#     plt.ylim(0., np.mean(np.log10(hist_log))*2.00)
#     plt.savefig(plot_direc+"Localization.IV.loglog."+filen+".%i.png" % indx, bbox_inches='tight')
#     plt.clf()
