#!/usr/bin/env python

import numpy as np
import re
from Scattering import scattering2D as scat


def read_xml(dir="./", file="Abmessung.xml"):
   '''
   Reads in 'param' objects from xml and
   return dictionary containing entries as
   strings.
   '''

   params = dict()

   f = open(dir + file, 'r')
   for line in f.readlines():
      words = filter(None, re.split('[ \n"<>]', line))
      if (words and words[0] == 'param'):
         params[words[2]] = words[3]

   return params


def write_states(states,
                 nr,
                 n_states,
                 nin_Max,
                 energ,
                 filen,
                 append=0,
                 state_sort=''):
   '''
   write states to file in proper format for Flo's code

   parameters:
   states: ndarray, float
      matrix containing row-wise states to write
   nr: int
      number of different energies = number of data sets
   n_states: int
      number of states to be propagated by Flo's code
   nin_Max: int
      actual number of open modes
   energ: ndarray, float
      scattering energy (k**2/2)
   filen: str
      namestring of calculation
   append: int
      whether coeffs should be appended to existing file
   state_sort: str (optional)
      string containing information which kind of
      states are calculated and to be propagated
   '''

   if (state_sort == ''):
      if (append):
         f = open('coeff.'+filen+'.dat', 'a')
         f.write('\n')
      elif(not append):
         f = open('coeff.'+filen+'.dat', 'w')
         f.write('%i -1\n' % nr)
         f.write('1.0 0.5\n')
   else:
      if (append):
         f = open('coeff.'+filen+'.'+state_sort+'.dat', 'a')
         f.write('\n')
      elif(not append):
         f = open('coeff.'+filen+'.'+state_sort+'.dat', 'w')
         f.write('%i -1\n' % nr)
         f.write('1.0 0.5\n')

   states = scat.inflate_small_mat(states, nin_Max)
   # append zeros to make number and size of vectors match
   # nin_Max = actual number of open modes in system


   f.write(str(energ)+' '+str(nin_Max)+'\n')
   for i in range(nin_Max):
      f.write('\n')
      if (i < n_states):
         f.write('1.0\n')
      else:
         f.write('2.0\n')
      for j in range(nin_Max):
         f.write('(' + str(states[i,j].real)
                 + ', ' + str(states[i,j].imag) + ')\n')
   return


def read_input(location, filen='pinput.dat'):
    '''
    Read in file 'pinput.dat', where data is stored in two columns:
    name of parameter, value of parameter
    returns dictionary of names and associated values, both as strings
    '''
    data = np.loadtxt(location + '/' + filen, dtype='string', unpack=True)
    try:
        params = dict(zip(data[0], data[1]))
        return params
    except:
        print "ABORT: error in reading input"
        return {}
  

def sort_eig(A):
    """
    Wrapper for np.linalg.eig(A) that returns:

    *) eigenvalues of A sorted in ascending order
    wrt their real values

    *) eigenvectors in the proper ordering (columnwise)
    """

    # get eigenvalues, eigenvectors of A
    eigenValues, eigenVectors = np.linalg.eig(A)  

    idx = np.absolute(eigenValues[:]).argsort()
    eigenValues[:] = eigenValues[idx]
    eigenVectors[:,:] = eigenVectors[:,idx]
    
    return eigenValues, eigenVectors

