import numpy as np

def calc_phase_times(S, dims, dk):
    '''
    Calculates phase times for all scattering
    amplitudes stored in S for a single radius
    and configuration.
    '''
    Nk = np.shape(dims)[0]
    pht = []
    for w in range(Nk):
        phi = np.absolute(S[3*w+1])**2 * (np.angle(S[3*w+2]) - np.angle(S[3*w+0])) / (2*dk)
        #phi = (np.angle(S[3*w+2]) - np.angle(S[3*w+0])) / (2*dk)
        pht.append(phi)

    return np.array(pht)


def calc_phase_means(pht, dims):
    '''
    Calculates mean phase times for a single
    radius and configuration.
    '''
    Nk = np.shape(pht)[0]
    means = np.zeros((Nk), dtype='complex')
    for w in range(Nk):
        means[w] = 1./(2.*dims[w]) * np.sum(pht[w])

    return means

