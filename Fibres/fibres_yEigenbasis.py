#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import pickle

from Scattering import scattering2D as scat
from Utils import utils as ut

par = ut.read_input('.')
try:
   filen = par['filen'] 
   lead_width = float(par['lead_width'])
   nyout = int(par['nyout'])
   modes_min = float(par['modes_min'])
   modes_max = float(par['modes_max'])
   nin = int(par['refmodes']) # nin = refmodes
   ptc = int(par['ptc'])
except KeyError:
    raw_input("WARNING: parameter missing in pinput.dat")

kmean = 0.5 * (modes_max + modes_min) * np.pi / lead_width
dk =  kmean * 10**(-8)
nin_Max = int(0.5 * (modes_max + modes_min)) # nin_Max = n_open_modes
dy = lead_width / (nyout + 1)
yrange = np.arange(1,nyout+1)*dy

S, dims, energs = scat.read_S('./', filen)

t = []
t_q = []
for i in range(3):
    m = dims[i]/2
    n = min(m, nin)
    t.append( S[i][m:,0:m] )      # full t-matrix
    t_q.append( S[i][m:m+n,0:n] ) # considered part of t-matrix 

def calc_modes(n):
   '''
   calculate normalized n-th lead mode

   In Flo's code y-axis points downwards, y=0 is upper egde of lead;
   modes are chosen such that most upper point is real and positive.
   '''
   chi_n = np.sin((n+1)*np.pi/lead_width * yrange)
   return np.sqrt(2./lead_width) * chi_n

chi = np.array(map(calc_modes, np.arange(0,nin)))

def calc_yoperator(chi):
   '''calculate operator y-d/2 in chi-representation'''
   return np.dot(chi.conj() * [list(yrange),]*nin, chi.T)*dy

y_op = calc_yoperator(chi)
y_eigm = ut.sort_eig(y_op)[1].T
#y_peaks = np.dot(y_eigm, chi)
y_peaks = np.dot(y_eigm, chi)



ut.write_states(y_eigm, nin, 'ypeaks')

q = scat.calc_Q(t_q, dk, inv=True)[0]
q_y = np.dot(y_eigm.conj(), np.dot(q, y_eigm.T))

t_y = np.dot(y_eigm.conj(), np.dot(t_q[1], y_eigm.T))

pos_range = (0, int(0.33*nin-0.5*ptc), int(0.5*nin-0.5*ptc))#, int(0.67*nin-0.5*ptc), nin-ptc)
# pos_range: place states to these positions

q_peaks = np.zeros((ptc,nin), dtype='complex')
for i, pos in enumerate(pos_range):
   peak_range = range(pos,pos+ptc)
   q_block = q_y[np.meshgrid(peak_range,peak_range)]

   block_eigv = ut.sort_eig(q_block)[1].T
   q_peaks = np.dot(block_eigv, y_eigm[peak_range])

   block_eigv_inv = np.linalg.inv(block_eigv).T
   y_eigm_inv = np.linalg.inv(y_eigm).T
   q_Lpeaks = np.dot(block_eigv_inv, y_eigm_inv[peak_range])

   pickle.dump((q_peaks, q_Lpeaks), open('qpeaks.%i.p'%i, 'wb'))
   ut.write_states(q_peaks, ptc, 'qpeaks.%i'%i)

tdt_peaks = np.zeros((ptc,nin), dtype='complex')
for i, posI in enumerate(pos_range):
   for j, posO in enumerate(pos_range):
      peak_rangeI = np.arange(posI,posI+ptc)
      peak_rangeO = np.arange(posO,posO+ptc)
      t_block = t_y[np.meshgrid(peak_rangeI,peak_rangeO)]
      #t_block = t_y[posO:posO+ptc, posI:posI+ptc]
      tdt_block = np.dot(t_block.T.conj(), t_block)

      block_eigv = ut.sort_eig(tdt_block)[1].T
      tdt_peaks = np.dot(block_eigv, y_eigm[peak_rangeI])

      block_eigv_inv = np.linalg.inv(block_eigv).T
      y_eigm_inv = np.linalg.inv(y_eigm).T
      tdt_Lpeaks = np.dot(block_eigv_inv, y_eigm_inv[peak_rangeI])

      pickle.dump((tdt_peaks, tdt_Lpeaks), open('tdtpeaks.%i.%i.p'%(i,j), 'wb'))
      ut.write_states(tdt_peaks, ptc, 'tdtpeaks.%i.%i'%(i,j))

#for n in range(nin):
   #plt.plot(yrange, chi[n])
#   plt.plot(yrange, y_peaks[n])
#   plt.show()
