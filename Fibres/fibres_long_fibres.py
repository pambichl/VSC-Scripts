#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
import sys
import os
from Utils import utils as ut
from Packets import transmission as trans
from Scattering import scattering2D as scat


Pi = np.pi
I = 0. + 1.J


parser = OptionParser()
parser.add_option("-f", "--filen",
                  action="store",
                  type = "str",
                  dest = "filen")
parser.add_option("-i", "--filemin",
                  action="store",
                  type = "int",
                  dest = "filemin")
parser.add_option("-a", "--filemax",
                  action="store",
                  type = "int",
                  dest = "filemax")
parser.add_option("-n", "--nseg",
                  action="store",
                  type = "int",
                  dest = "Nseg")
parser.add_option("-p", "--npot",
                  action="store",
                  type = "int",
                  dest = "Npot")

options = parser.parse_args()[0]
filen = options.filen
filemin = options.filemin
filemax = options.filemax
Nseg = options.Nseg
Npot = options.Npot

if (not filen): sys.exit("ABORT: specify filename (str, -f, --filen)")
if (not filemin): sys.exit("ABORT: specify minimum calculation identifier of S-matrices (int, -i, --filemin)")
if (not filemax): sys.exit("ABORT: specify maximum calculation identifier of S-matrices (int, -a, --filemax)")
if (not Nseg and Nseg!=0):
   print "WARNING: number of segments in each block not specified (int, -n, --nseg)"
   Nseg = 5
if (not Npot and Npot!=0):
   print "WARNING: number of potentiating steps not specified (int, -p, --npot)"
   Npot = 4

Nfile = filemax - filemin + 1

data_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
              "VSC2/Fibres-20130802/" + filen + "/")

Slist = []

for i in np.arange(filemin, filemax+1):

   file_i = "fi%i" % i
   scatter_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                 "VSC2/Fibres-20130802/fi%i/" % i)

   print "Reading scattering matrices fi%i" % i
   if (not os.path.exists("scatter."+file_i+".dat.npy")):
      S, dims, energs = scat.read_S(scatter_direc, file_i, old_ver=0)
      np.save(data_direc + "scatter." + file_i + ".dat", (S, dims, energs))
   else:
      S, dims, energs = np.load(data_direc + "scatter." + file_i + ".dat.npy")
   
   Slist.append(S)


for i in range(max(0,Nseg-Nfile)):
   el = np.random.randint(low=0, high=len(Slist))
   Slist.append(Slist[el]) # if Nseg > Nfile append randomly chosen S-matrices


Slist = np.array(Slist)
Nw = len(dims)
dw = np.sqrt(2.0*energs[2]) - np.sqrt(2.0*energs[1])

tsafe = trans.calc_t(Slist[0], dims, 1)[0]


def compose_block(Slist, dims, Nseg):
   '''
   Compose block of size Nseg from segments stored in Slist.
   '''
   
   comp = np.empty((Nseg,), dtype="int")
   if (Nseg > 1):
      seg_ind = np.random.randint(low=0, high=Nseg)
      Seg = Slist[seg_ind] # choose a segment from list as starting point
      comp[0] = seg_ind+1
   else:
      Seg = Slist[0]
      comp[0] = 1
      
   for seg in range(Nseg-1):

      t, _ = trans.calc_t(Seg, dims, 1)
      tp, _ = trans.calc_tp(Seg, dims, 1)
      r, _ = trans.calc_r(Seg, dims, 1)
      rp, _ = trans.calc_rp(Seg, dims, 1)
         
      flip = np.random.choice([True,False]) # flip segment or not

      if (Nseg > 1):
         seg_ind = np.random.randint(low=0, high=Nseg)
         SA = Slist[seg_ind] # choose a segment to add
         if(not flip): comp[seg+1] = seg_ind+1
         if(flip): comp[seg+1] = -(seg_ind+1)
      else:
         SA = Slist[0]
         if(not flip): comp[seg+1] = 1
         if(flip): comp[seg+1] = -1
         
      if (not flip):
         tA, _ = trans.calc_t(SA, dims, 1)
         tpA, _ = trans.calc_tp(SA, dims, 1)
         rA, _ = trans.calc_r(SA, dims, 1)
         rpA, _ = trans.calc_rp(SA, dims, 1)
      elif (flip):
         tA, _ = trans.calc_tp(SA, dims, 1)
         tpA, _ = trans.calc_t(SA, dims, 1)
         rA, _ = trans.calc_rp(SA, dims, 1)
         rpA, _ = trans.calc_r(SA, dims, 1)

      Seg = [] # new S matrix for old + new at each frequency
      Nw = len(dims)
      for w in range(Nw):

         tN = np.dot(tA[w], np.dot(np.linalg.inv(np.eye(dims[w]/2) - np.dot(rp[w], rA[w])), t[w]))
         tpN = np.dot(tp[w], np.dot(np.linalg.inv(np.eye(dims[w]/2) - np.dot(rA[w], rp[w])), tpA[w]))
         rN = r[w] + np.dot(tp[w], np.dot(np.linalg.inv(np.eye(dims[w]/2) - np.dot(rA[w], rp[w])), np.dot(rA[w], t[w])))
         rpN =  rpA[w] + np.dot(tA[w], np.dot(np.linalg.inv(np.eye(dims[w]/2) - np.dot(rp[w], rA[w])), np.dot(rp[w], tpA[w])))
         
         s = np.zeros((dims[w], dims[w]), dtype="complex")
         s[:dims[w]/2,:dims[w]/2] = rN
         s[:dims[w]/2,dims[w]/2:] = tpN
         s[dims[w]/2:,:dims[w]/2] = tN
         s[dims[w]/2:,dims[w]/2:] = rpN
      
         Seg.append(s)

   return Seg, comp


def create_output(S, energs, filen):
   '''
   Creates file Smat.filen.dat.
   '''

   Nw = len(energs)
   with open("Smat."+filen+".dat", "w") as f_handle:
      f_handle.write('%i\n\n'%Nw)

      for w in range(Nw):
         dim = np.shape(S[w])[0]
         f_handle.write('%.12f\n'%energs[w])
         f_handle.write('%i\n\n'%dim)
         for i in range(dim):
            for j in range(dim):
               f_handle.write('%3i %3i  %20.16f  %20.16f \n' % (i, j, S[w][i,j].real, S[w][i,j].imag))
            f_handle.write('\n')


for pot in range(Npot):

   Seglist = [] # segments blocks will be randomly built of
   print "Calculating step %i of %i" % (pot+1,Npot)
   
   if (pot==Npot-1):
      Nblock=1
   else:
      Nblock = Nseg

   for block in range(Nblock):

      print "Composing block %i of %i" % (block+1,Nseg),
      Seg, comp = compose_block(Slist, dims, Nseg)
      Seglist.append(Seg) # get new segments by composing old segments to blocks
      print comp

   Slist = Seglist # define blocks as new segments

Slist = np.array(Slist)

S = Slist[0]


print
print "Transmission eigenvalues at center frequency:"
print np.abs(ut.sort_eig(trans.calc_t(S, dims, 1)[0][Nw/2])[0])
print "Reflection eigenvalues at center frequency: "
print np.abs(ut.sort_eig(trans.calc_r(S, dims, 1)[0][Nw/2])[0])

Q = scat.calc_Q((S[3*(Nw/3/2)+0],S[3*(Nw/3/2)+1],S[3*(Nw/3/2)+2]), dw, inv=False)[0]
Q = scat.calc_Q((S[0],S[1],S[2]), dw, inv=False)[0]
print "Time delay eigenvalues at center frequency:"
print ut.sort_eig(Q)[0]

print "Writing output"
create_output(S, energs, filen)





