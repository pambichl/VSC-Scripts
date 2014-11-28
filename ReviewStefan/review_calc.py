#!/usr/bin/env python

import numpy as np

from Scattering import scattering2D as scat
from Utils import utils as ut
import sys

if (len(sys.argv)!=2):
    sys.exit("ABORT: filen needed")

filen = str(sys.argv[1]) # namestring of calculation

S_dis, dims, energs = scat.read_S('./', filen+'.dis')
S_clean, dims, energs = scat.read_S('./', filen+'.clean')

dk = np.sqrt(2.0*energs[-1]) - np.sqrt(2.0*energs[-2])

t_dis = []
t_clean = []
for i in range(3):
    nin = dims[i]/2
    t_dis.append( S_dis[i][nin:,0:nin] )
    t_clean.append( S_clean[i][nin:,0:nin] )

Q = scat.calc_Q(S_clean,dk)[0]
Q11 = Q[0:nin,0:nin]
Q11EigVec = ut.sort_eig(Q11)[1]

TransStates = np.dot( np.linalg.inv(t_dis[1]), np.dot( t_clean[1], Q11EigVec ))

ut.write_states(Q11EigVec.T, nin, nin, energs[-2], filen+'.clean')
ut.write_states(TransStates.T, nin, nin, energs[-2], filen+'.dis')

