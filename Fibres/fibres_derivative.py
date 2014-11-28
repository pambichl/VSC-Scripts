#!/usr/bin/env python

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from cmath import *
import matplotlib.pyplot as plt
import sys

import pickle

from Utils import utils as ut
from Scattering import scattering2D as scat

from Packets import transmission as trans

Pi = np.pi
I = 0. + 1.J

par = ut.read_input('.')
try:
   filen = par['filen'] # namestring of calculation where scatter data is stored
   pic_filen = par['pic_filen'] # namestring of calculation where waveplots are stored
   lead_width = float(par['lead_width'])
   modes_min = float(par['modes_min'])
   modes_max = float(par['modes_max'])
   refmodes = int(par['refmodes']) # number of modes to be considered
except KeyError:
    raw_input("WARNING: parameter missing in pinput.dat")

expo = 1

kmean = 0.5 * (modes_max + modes_min) * np.pi / lead_width
dk =  kmean * 10**(-expo)
nin_Max = int(0.5 * (modes_max + modes_min)) # nin_Max = n_open_modes



try:
    # on which machine calculation was performed, VSC or VSC2
    machine = int(sys.argv[1])
except IndexError:
    machine = 2

if machine == 1:
    data_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                  "VSC/Fibres-20130819/" + filen + "/")
    pic_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                 "VSC/Fibres-20130819/" + pic_filen + "/")
elif machine ==2:
    data_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                  "VSC2/Fibres-20130802/" + filen + "/")
    pic_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                 "VSC2/Fibres-20130802/" + pic_filen + "/")
else:
    print "specify machine (VSC: 1, VSC2: 2)"
plot_direc = data_direc


S, dims = scat.read_S(data_direc, filen+".1E-%i"%expo, old_ver=1)
modes = int(max(dims)/2)


t , t_q = trans.calc_t(S, dims, refmodes)

tp = (t[2] - t[0]) / (2*dk)

ShiftPhase = -I * (t[1]/np.absolute(t[1])).conj() * (t[2]/np.absolute(t[2]) - t[0]/np.absolute(t[0])) / (2*dk)

#ShiftPhase = 0.001*np.ones((3,3), dtype="complex")

ShiftFacsP = np.exp( I * ShiftPhase * dk)
ShiftFacsM = np.exp(-I * ShiftPhase * dk)

#print np.absolute(ShiftFacsM)

tpShift = (t[2]*ShiftFacsM - t[0]*ShiftFacsP) / (2*dk) + I * t[1]*ShiftPhase

# print np.absolute(t[1]) * (t[2]/np.absolute(t[2]) - t[0]/np.absolute(t[0])) / (2*dk)\
#     + (np.absolute(t[2]) - np.absolute(t[0])) / (2*dk) * t[1]/np.absolute(t[1]) # Kettenregel
# print tp # normal
# print tpShift # mit shift

q = -I * np.linalg.inv(t[1]) * tp
print np.linalg.eig(q)[0]
qShift = -I * np.linalg.inv(t[1]) * tpShift 
print np.linalg.eig(qShift)[0]

# q0 = scat.calc_Q(t, dk, inv=True)[0]
# shift = np.array([np.exp(I*q0*dk), np.ones((3,3), dtype='complex'), np.exp(-I*q0*dk)])

# tS = shift * t
# qS = scat.calc_Q(tS, dk, inv=True)[0]

# #print t[1]
# #print tS[1]

# #print q0
# #print qS+q0

# print np.linalg.eig(q0)[0]
# print np.linalg.eig(qS-q0)[0]

# epsilon, V = np.linalg.eig(t)
# phi = epsilon / np.absolute(epsilon)
# q0 = -I * phi[1].conj() * (phi[2] - phi[0]) / (2*dk)

# # VL = np.linalg.inv(V)

# # tp = np.dot(V[1], np.dot(np.diag(epsilon[1]), VL[1]))

# # qp = -I * np.dot(V[1], np.dot( np.diag( 1.0/epsilon[1] * (epsilon[2] - epsilon[0]) / (2*dk) ), VL[1]))
# # #qp = -I * np.dot(V[1], np.dot( np.diag( phi[1].conj() * (phi[2] - phi[0]) / (2*dk) ), VL[1]))


# # q = scat.calc_Q(t, dk, inv=True)[0]
  
# # #print tp
# # #print t[1]

# # print qp
# # print q

# # print np.linalg.eig(qp)[0]
# # print np.linalg.eig(q)[0]

# VL = np.linalg.inv(V)
# epsilon_mat = np.array(map(np.diag, epsilon))
# #q0_mat = np.eye(3, dtype="complex")
# pref = [1.0, 0.0, -1.0]
# q0_mat = np.array([np.diag(np.exp(pref[i] * I * q0 * dk)) for i in range(3)])

# #print q0_mat

# tShift = np.array([np.dot(V[i], np.dot(epsilon_mat[i], np.dot(q0_mat[i], VL[i]))) for i in range(3)])
# #print tShift - np.array(t)

# q = scat.calc_Q(t, dk, inv=True)[0]
# qShift = scat.calc_Q(tShift, dk, inv=True)[0]

# print np.linalg.eig(q)[0]
# print q0.real
# print (np.linalg.eig(qShift)[0] + q0).real

