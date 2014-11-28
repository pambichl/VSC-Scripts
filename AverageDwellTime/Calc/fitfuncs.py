import numpy as np


def gauss(x, x0, s, a):
    return a * np.exp(-(x-x0)**2/(2.*s))

def diffusive(T, a):
    '''transmission distribution for diffusive system'''
    T0 = 0.01
    return a / (T * np.sqrt(1.-T)) * 0.5 * (np.sign(T-T0)+1.)

def chaotic(T, a):
    '''transmission distribution for chaotic system'''
    return a / np.sqrt(T * (1.-T))

def localized_log_amp_only(T, mean_lnG):
    '''total transmission distribution for 1D localized system'''
    return (T * np.sqrt(np.arccosh(T**(-1./2.))) / (T**(3./2.) * (1.0-T)**(1./4.)) * np.exp(1./mean_lnG * np.arccosh(T**(-1./2.))**2)).real

def localized(T, a, xi):
    '''
    total transmission distribution for 1D localized system;
    effective localization length also fit parameter (in units
    of L)
    '''
    return (a * np.sqrt(np.arccosh(T**(-1./2.))) / (T**(3./2.) * (1.-T)**(1./4.)) \
        * np.exp(-0.5 * xi * np.arccosh(T**(-1./2.))**2.)).real

def localized_log(T, a, xi):
    '''
    total transmission distribution for 1D localized system;
    distribution for lnT ( P(lnT) = T * P(T); prefactor ln(10), which can
    be considered as contained in amplitude a
    effective localization length also fit parameter (in units
    of L)
    '''
    return (a * T * np.sqrt(np.arccosh(T**(-1./2.))) / (T**(3./2.) * (1.-T)**(1./4.)) \
        * np.exp(-0.5 * xi * np.arccosh(T**(-1./2.))**2.)).real
