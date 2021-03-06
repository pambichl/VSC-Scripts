#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import pickle

from Scattering import scattering2D as scat
from Utils import utils as ut

par = ut.read_input('.')
try:
   filen = par['filen'] 
   modes_min = float(par['modes_min'])
   modes_max = float(par['modes_max'])
   lead_width = float(par['lead_width'])
   nin = int(par['refmodes']) # nin = refmodes
   Nk = int(par['Nk'])
except KeyError:
    raw_input("WARNING: parameter missing in pinput.dat")

kmean = 0.5 * (modes_max + modes_min) * np.pi / lead_width
dk =  kmean * 10**(-8)
Dk = (modes_max - modes_min) * np.pi / lead_width / (Nk -1)

S, dims, energs = scat.read_S('./', filen)
if (len(energs)==0):
   energs = np.linspace(modes_min, modes_max, Nk) * np.pi/lead_width
   energs = 0.5 * energs**2
else:
   energs = energs[::3]

centpos = int(0.5*Nk)
quartpos = int(0.25*Nk)
nin_Max = int(0.5*np.max(dims[::3][centpos])) # nin_Max = n_open_modes

t = []
t_q = []
for i in np.arange(3*Nk):
   m = dims[i]/2
   n = min(m, nin)
   t.append( S[i][m:,0:m] )      # full t-matrix
   t_q.append( S[i][m:m+n,0:n] ) # considered part of t-matrix 

q = []
for i in range(Nk):
    q.append(scat.calc_Q([t_q[3*i+0],t_q[3*i+1],t_q[3*i+2]], dk, inv=True)[0])
q = np.array(q)

t_q = t_q[::3]
pos = ((-1,0), (-quartpos-1,quartpos), (-centpos,centpos-1))
pos_string = ('full', 'quart', 'cent')

for i, p, s in zip(range(3), pos, pos_string):
   KEigVal, KEigVec = ut.sort_eig( np.dot(np.linalg.inv(t_q[p[1]]), t_q[p[0]]) )[:2]
   KEigVec = KEigVec.T
   KEigVec_L = np.linalg.inv(KEigVec).T
   data = np.array([range(nin), (-1.0J*np.log(KEigVal)).real, (-1.0J*np.log(KEigVal)).imag])
   np.savetxt('KEigVals.'+filen+'.'+s+'.txt', data.T, fmt=('%3i', '%12.8f', '%12.8f'))
   ut.write_states(KEigVec, nin, nin_Max, energs[centpos], filen, s)
   pickle.dump((KEigVec, KEigVec_L), open('KEigVecs.%i.p'%i, 'wb'))




