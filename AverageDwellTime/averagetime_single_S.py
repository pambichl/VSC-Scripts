#!/usr/bin/env python

import numpy as np
import scipy as sp
import scipy.optimize as spopt
import matplotlib.pyplot as plt
import sys
from math import pi as Pi

from Utils import utils as ut
from Scattering import scattering2D as scat

filen = str(sys.argv[1]) # namestring of calculation
n     = int(sys.argv[2]) # number of calculation
dims  = int(sys.argv[3]) # number of open modes in left lead

data_direc = "../../VSC2/AverageDwellTime-20130618/" + filen + "/scatterdata/"

energs = scat.read_S(data_direc, filen+".%i.0"%n)[2]
kvals = np.sqrt(2*energs)
dk = (kvals[-1] - kvals[-2]) / 1.0

filen_i = filen+".%i.%i" % (n,0)
print "reading file: "+filen_i
S = scat.read_S(data_direc, filen_i)[0]
q_mean = np.trace(scat.calc_Q(S,dk)[0]) /(2*dims)

print q_mean.real
