#!/usr/bin/env python

import numpy as np
from Scattering import scattering2D as scat
from Utils_green import utils as ut
import matplotlib.pyplot as plt

filen = "sp10"
file_path = "/home/lv70071/ambichl/greens_code/calculations/Spiral-20141105/" + filen + "/"


S, dims, energs = scat.read_S(file_path, filen, old_ver=0)
nin = dims/2
Nw = np.shape(dims)[0]

t = []
tdt = []
for w in range(Nw):
    t.append(S[w][:nin[w],nin[w]:])
    tdt.append(t[-1].conj().T.dot(t[-1]))

Nin = np.amax(nin)
tdt_EigVal = np.zeros((Nw,Nin), dtype="complex")
tdt_EigVec = np.zeros((Nw,Nin,Nin), dtype="complex")

for w in np.arange(0,Nw):
    tdt_EigVal[w][:nin[w]] = ut.sort_eig(tdt[w])[0]
    tdt_EigVec[w][:nin[w],:nin[w]] = ut.sort_eig(tdt[w])[1].T
    states = tdt_EigVec[w]
    if (w == 0):
       ut.write_states(states, Nw, Nin, Nin, energs[0], filen, append=0)
    ut.write_states(states, Nw, Nin, Nin, energs[w], filen, append=1)

def calc_int_corr(vecs, i0):
    Nw = np.shape(vecs)[0]
    ints = np.absolute(vecs)**2
    corrs = np.zeros((Nw,), dtype="float")
    for w in range(Nw):
        corrs[w] = np.mean(ints[i0] * ints[w]) / (np.mean(ints[i0]) * np.mean(ints[w])) - 1 

    return corrs / corrs[i0]

def calc_prop(t, vec):
    Nw = np.shape(t)[0]
    out = np.zeros((Nw, Nin), dtype="complex")
    for w in range(Nw):
        out[w,:nin[w]] = t[w].dot(vec[:nin[w]])
    return out

def calc_modes(k):
    yvals = np.arange(0.0, 1.0, 0.01)
    ky = np.arange(1,Nin+1,1) * np.pi / 1.0
    kx = np.sqrt(np.abs(k**2 - ky**2))
    flux = 1.0 / np.sqrt(kx)
    chi = np.sqrt(2.0 / 1.0) * np.sin(np.einsum('y,n->yn', yvals, ky))
    return  np.einsum('n,yn->yn', flux, chi)

yvals = np.arange(0.0, 1.0, 0.01)
### state to inspect
i0 = 7
outvec = calc_prop(t, tdt_EigVec[i0][np.argmax(tdt_EigVal[i0])])
###
real_space_out = np.zeros((Nw, len(yvals)), dtype="complex")
for w in range(Nw):
    k = np.sqrt(2.0 * energs[w])
    real_space_out[w] = np.einsum('n,yn->y', outvec[w], calc_modes(k))
    
fig = plt.figure(0)
plt.ylabel("frequency")
plt.xlabel("spatial transmitted intensity")
plt.pcolor(np.reshape(np.absolute(real_space_out)**2, (Nw, len(yvals))))

plt.figure(1)
plt.title("reference frequency index i0=%i" % i0)
plt.xlabel("frequency")
plt.ylabel("transmission")
plt.plot(map(np.linalg.norm, outvec))

plt.figure(2)
plt.plot(calc_int_corr(outvec, i0), 'b')
plt.plot(calc_int_corr(outvec, i0), 'g')
plt.plot(calc_int_corr(outvec, i0), 'r')

plt.show()

