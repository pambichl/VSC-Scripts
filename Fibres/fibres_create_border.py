#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
from Utils import utils as ut

parser = OptionParser()
parser.add_option("-d", "--dir",
                  action="store",
                  type = "str",
                  dest = "dir")
parser.add_option("-f", "--file",
                  action="store",
                  type = "str",
                  dest = "file")


options = parser.parse_args()[0]
dir = options.dir
file = options.file

if(not dir): dir="./"
if(not file): file="Abmessung.xml"

params = ut.read_xml(dir, file)

modes_min = float(params['modes_min'])
modes_max = float(params['modes_max'])
pphwl = int(params['points_per_hwl'])
W = float(params['W'])
L = float(params['L'])
amp = float(params['borderamp'])

ny = 0.5 * (modes_min + modes_max) * pphwl
dx = W / (ny+1.0)
nx = L / dx

print "Writing border."

x=np.linspace(0.,L,nx)
f=amp*np.cos(x*np.pi/L)
np.savetxt('border.dat',zip(x,f))
