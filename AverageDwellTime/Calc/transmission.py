import numpy as np
import scipy as sp
from Utils import utils as ut

def calc_t(S, dims, (nr,nfreq,nconf,nin_max), fullout=False):
    '''calculate t-matrices from S-matrices'''
    if (not fullout):
        t = np.zeros((nr, nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(nfreq):
                for k in range(nconf):
                    t_mat = np.array(S[i,k,1+3*j][dims[j]:2*dims[j],0:dims[j]])
                    # if transmission so low that somewhere set to zero:
                    zero_mask = t_mat == 0.
                    t_mat[zero_mask] = 10**(-16)
                    t[i,j,k,0:dims[j],0:dims[j]] = t_mat
    elif (fullout):
        t = np.zeros((nr, 3*nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(3*nfreq):
                for k in range(nconf):
                    t_mat = np.array(S[i,k,j][dims[int(j/3)]:2*dims[int(j/3)],0:dims[int(j/3)]])
                    # if transmission so low that somewhere set to zero:
                    zero_mask = t_mat == 0.
                    t_mat[zero_mask] = 10**(-16)
                    t[i,j,k,0:dims[int(j/3)],0:dims[int(j/3)]] = t_mat
    return t


def calc_r(S, dims, (nr,nfreq,nconf,nin_max), fullout=False):
    '''calculate r-matrices from S-matrices'''
    if (not fullout):
        r = np.zeros((nr, nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(nfreq):
                for k in range(nconf):
                    r_mat = np.array(S[i,k,1+3*j][0:dims[j],0:dims[j]])
                    # if reflection so low that somewhere set to zero:
                    zero_mask = r_mat == 0.
                    r_mat[zero_mask] = 10**(-16)
                    r[i,j,k,0:dims[j],0:dims[j]] = r_mat
    elif (fullout):
        r = np.zeros((nr, 3*nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(3*nfreq):
                for k in range(nconf):
                    r_mat = np.array(S[i,k,j][0:dims[int(j/3)],0:dims[int(j/3)]])
                    # if reflection so low that somewhere set to zero:
                    zero_mask = r_mat == 0.
                    r_mat[zero_mask] = 10**(-16)
                    r[i,j,k,0:dims[int(j/3)],0:dims[int(j/3)]] = r_mat
    return r


def calc_tp(S, dims, (nr,nfreq,nconf,nin_max), fullout=False):
    '''calculate t-matrices from S-matrices'''
    if (not fullout):
        tp = np.zeros((nr, nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(nfreq):
                for k in range(nconf):
                    tp_mat = np.array(S[i,k,1+3*j][0:dims[j],dims[j]:2*dims[j]])
                    # if transmission so low that somewhere set to zero:
                    zero_mask = tp_mat == 0.
                    tp_mat[zero_mask] = 10**(-16)
                    tp[i,j,k,0:dims[j],0:dims[j]] = tp_mat
    elif (fullout):
        tp = np.zeros((nr, 3*nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(3*nfreq):
                for k in range(nconf):
                    tp_mat = np.array(S[i,k,j][0:dims[int(j/3)],dims[int(j/3)]:2*dims[int(j/3)]])
                    # if transmission so low that somewhere set to zero:
                    zero_mask = tp_mat == 0.
                    tp_mat[zero_mask] = 10**(-16)
                    tp[i,j,k,0:dims[int(j/3)],0:dims[int(j/3)]] = tp_mat
    return tp


def calc_rp(S, dims, (nr,nfreq,nconf,nin_max), fullout=False):
    '''calculate t-matrices from S-matrices'''
    if (not fullout):
        rp = np.zeros((nr, nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(nfreq):
                for k in range(nconf):
                    rp_mat = np.array(S[i,k,1+3*j][dims[j]:2*dims[j],dims[j]:2*dims[j]])
                    # if transmission so low that somewhere set to zero:
                    zero_mask = rp_mat == 0.
                    rp_mat[zero_mask] = 10**(-16)
                    rp[i,j,k,0:dims[j],0:dims[j]] = rp_mat
    elif (fullout):
        rp = np.zeros((nr, 3*nfreq, nconf, nin_max, nin_max), dtype='complex')
        for i in range(nr):
            for j in range(3*nfreq):
                for k in range(nconf):
                    rp_mat = np.array(S[i,k,j][dims[int(j/3)]:2*dims[int(j/3)],dims[int(j/3)]:2*dims[int(j/3)]])
                    # if transmission so low that somewhere set to zero:
                    zero_mask = rp_mat == 0.
                    rp_mat[zero_mask] = 10**(-16)
                    rp[i,j,k,0:dims[int(j/3)],0:dims[int(j/3)]] = rp_mat
    return rp


def calc_tdt_eigvals(t, (nr,nfreq,nconf,nin_max)):
    '''Calculate eigenvalues of t^dagger.t from t-matrices'''
    tdt_eigvals = np.zeros((nr, nfreq, nconf, nin_max), dtype='complex')
    for i in range(nr):
        for j in range(nfreq):
            for k in range(nconf):
                tdt_eigvals[i,j,k] = ut.sort_eig(np.dot(t[i,j,k].T.conj(), t[i,j,k]))[0]
    # sometimes transmission so low that somewhere tdt-eigval gets negative
    neg_mask = tdt_eigvals < 10**(-16)
    tdt_eigvals[neg_mask] = 10**(-16)
    return tdt_eigvals
