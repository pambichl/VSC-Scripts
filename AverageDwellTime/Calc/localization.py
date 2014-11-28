import numpy as np


def check_single_channel(tdt_eigvals, n, (nr,nfreq,nconf)):
    '''
    Print highest transmission eigenvalue and ratios of second and third highest to
    highest transmission eigenvalue for first n configurations at center frequency.
    If single-channel transport, ratios should be very small.

    int n: number of configurations to print
    '''
    j = int(nfreq/2)
    for i in range(nr):
        print "highest transmission and ratios of transmissions at radius #%i:"%i
        for k in range(min(n,nconf)):
            eigv = np.sort(tdt_eigvals[i,j,k]).real
            print "%5f, %5f, %5f" % (eigv[-1], eigv[-2]/eigv[-1], eigv[-3]/eigv[-1])

def calc_localization_meas_I(t, dims, (nr,nfreq,nin_max)):
    '''
    Caculate measure for localization as used in Chabanov PRL 87
    using variance of t-matrix elements. If measure > 7/3, wave
    should be localized.
    '''
    T = np.abs(t)**2
    s_mean = np.zeros((nr, nfreq, nin_max, nin_max), dtype='float')
    for i in range(nr):
        for j in range(nfreq):
            for m in range(dims[j]):
                for n in range(dims[j]):
                    s_mean[i,j,m,n] = np.var(T[i,j,:,m,n].flatten()) / np.mean(T[i,j,:,m,n].flatten())**2

    s_mmean = np.empty((nr, nfreq), dtype='float')
    for i in range(nr):
        for j in range(nfreq):
            s_mmean[i,j] = np.sum(s_mean[i,j,:,:].flatten()) / dims[j]**2

    return s_mmean

def calc_localization_meas_II(tdt_eigvals, nr, bins):
    '''
    Calculate quantity used for localization-measure from Beenakker review Eq. (156).
    T_n = 1./cosh(x_n)^2 -> x_n = arccosh( sqrt( 1./T_n ) )
    Returns distribution of x_n's over frequencies and configurations; if crystallized, wave
    should be localized.
    '''
    zero_mask = np.array(tdt_eigvals!=0.+0.j)
    loc_var = np.array([np.arccosh(np.sqrt(1./tdt_eigvals[i][zero_mask[i]])) for i in range(nr)])
    cryst_hists = np.array([np.histogram(loc_var[i], bins=bins, density=False)[0] for i in range(nr)])
    return loc_var, cryst_hists

def calc_localization_meas_III(tdt_eigvals, (nr,nfreq,nconf), bins, ln_cond_range):
    '''
    Calculate logarithm of total transmissions as measure for localization.
    Returns distribution over frequencies and configurations; if Gaussian, wave
    should be localized.
    '''
    ln_cond = np.log(np.sum(tdt_eigvals, axis=3)).reshape((nr, nfreq*nconf))
    lognorm_hists = np.array([np.histogram(ln_cond[i], bins=bins, range=ln_cond_range, density=False)[0] for i in range(nr)])
    return ln_cond, lognorm_hists

def calc_localization_meas_IV(tdt_eigvals, bins, tdt_range, log=False):
    '''
    Calculate distribution of t^dagger.t - eigenvalues to see whether bimodal or not.
    '''
    tdt_hists = []
    if (log==True):
        vals_range = (np.log10(tdt_range[0]), np.log10(tdt_range[-1]))
    else:
        vals_range = tdt_range

    for vals in tdt_eigvals:
        tdt_hists.append(np.histogram(vals.real, bins=bins, range=vals_range, density=True)[0])
    return np.array(tdt_hists)
