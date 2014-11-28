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

nin = 50
nloop = 10000

#np.random.seed(1786)

A = np.random.random(size=(nin,nin)) + I * np.random.random(size=(nin,nin))
A = 0.5 * (A + A.conj().T)
A = ut.sort_eig(A)[1]
A = np.dot( A.conj().T, np.dot( 100*np.diagflat(np.random.random(size=(nin))), A)) 
#A = A / np.linalg.norm(A)
C = np.random.random(size=(nin,nin)) + I * np.random.random(size=(nin,nin))
C = 0.5 * (C + C.conj().T)
C = ut.sort_eig(C)[1]
C = np.dot( C.conj().T, np.dot( np.diagflat(np.random.random(size=(nin))), C)) 
#C = C / np.linalg.norm(C)
D = np.random.random(size=(nin,nin)) + I * np.random.random(size=(nin,nin))
D = 0.5 * (D + D.conj().T)
D = ut.sort_eig(D)[1]
D = np.dot( D.conj().T, np.dot( np.diagflat(np.random.random(size=(nin))), D)) 
#D = D / np.linalg.norm(D)

phi = np.random.random(size=(nin)) + I * np.random.random(size=(nin))
phi = phi / np.linalg.norm(phi)


def calc_coeffs(phi, (A, C, D)):
   '''
   Calculates coefficients for eigenvalue
   equation, i.e., the expectation values
   of the matrices A, C, D with vector phi
   '''
   a = np.dot( phi.conj(), np.dot( A, phi))
   c = np.dot( phi.conj(), np.dot( C, phi))
   d = np.dot( phi.conj(), np.dot( D, phi))

   return a, c, d

def calc_goal_matrix((a, c, d), (A, C, D)):
   '''
   Calculates matrix of which the eigenstates
   are determined in each iteration step.
   '''
   G = -((a - c**2) / d**2) * D + (1.0 / d) * A - (2.0 * c / d) * C #- 0.*np.eye(nin, dtype='complex')

   return G


for i in range(nloop):

   a, c, d = calc_coeffs(phi, (A,C,D))
   G = calc_goal_matrix((a,c,d), (A,C,D))

   #phi = (ut.sort_eig(G)[1].T)[-1]
   phi = np.dot(G, phi)
   phi = phi / np.linalg.norm(phi)
   # Fixpunktgleichung mit lambda = 0 einfach angenommen
   # konvergiert schneller gegen selben Vektor wie
   # Eigenwertgleichung mit groesztem Eigenwert folgend
   # manchmal...


a, c, d = calc_coeffs(phi, (A,C,D))
G = calc_goal_matrix((a,c,d), (A,C,D))

print np.dot(G, phi) / phi
print np.dot(G, np.dot(G, phi)) / phi
#print ut.sort_eig(G)[0]

print 1.0 / d * (a - c**2)

phi = np.random.random(size=(nin)) + I * np.random.random(size=(nin))
phi = phi / np.linalg.norm(phi)

a, c, d = calc_coeffs(phi, (A,C,D))
G = calc_goal_matrix((a,c,d), (A,C,D))

print 1.0 / d * (a - c**2)

print np.linalg.det(G)
#print ut.sort_eig(G)[0]

