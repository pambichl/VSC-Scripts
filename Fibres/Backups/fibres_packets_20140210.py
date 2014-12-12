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

kmean = 0.5 * (modes_max + modes_min) * np.pi / lead_width
dk =  kmean * 10**(-8)
nin_Max = int(0.5 * (modes_max + modes_min)) # nin_Max = n_open_modes



try:
    # on which machine calculation was performed, VSC or VSC2
    machine = int(sys.argv[1])
except IndexError:
    machine = 2

if machine == 1:
    data_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                  "VSC/Fibres-20130819/" + filen + "/")
elif machine ==2:
    data_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                  "VSC2/Fibres-20130802/" + filen + "/")
else:
    print "specify machine (VSC: 1, VSC2: 2)"

if machine == 1:
    pic_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                  "VSC/Fibres-20130819/" + pic_filen + "/")
elif machine == 2:
    pic_direc = ("/home/ambichl/Universitaet/Scattering-Data/"
                  "VSC2/Fibres-20130802/" + pic_filen + "/")

plot_direc = data_direc

fullout = False


S, dims = scat.read_S(data_direc, filen, old_ver=1)
modes = int(max(dims)/2)

if (len(dims)%3 == 0):
    Nk = len(dims)/3
else:
    sys.exit("ABORT: # of S matrices not divisible by 3")

if (refmodes > int(min(dims)/2)):
    print ("WARNING: # of modes under consideration (%i)"
           "> # minimum number of modes (%i)") % (refmodes,int(min(dims)/2))


t , t_q = trans.calc_t(S, dims, refmodes)
#trans.write_teigvals(t_q, refmodes)
q = []
for i in range(Nk):
    q.append(scat.calc_Q([t_q[3*i+0],t_q[3*i+1],t_q[3*i+2]], dk, inv=True)[0])
q = np.array(q)

refpos = int(0.5*(Nk-1)) # k[refpos] = kmean for odd Nk

qref = q[refpos]
tref = t_q[3*refpos+1]
tinv = np.linalg.inv(tref)


class TRC:
   '''
   '''
   modes = 0
   refmodes = 0
   refpos = 0

   Nk = 0

   t = 0
   tref = 0
   tinv = 0

   q = 0
   qeigbas = 0
   qLeigbas = 0

   init = False

   @classmethod
   def class_init(cls, modes, refmodes, refpos,
                  Nk, t, tref, q):
      TRC.modes = modes
      TRC.refmodes = refmodes
      TRC.refmodes = refpos
      TRC.Nk = Nk
      TRC.t = t
      TRC.tref = tref
      TRC.tinv = np.linalg.inv(tref)
      TRC.q = q
      TRC.init = True

      qeigbas, qLeigbas = trans.calc_qeigbas(q, Nk, modes, refmodes)

   def __init__(self, op_char='m', label=r'modes', color='b',
                operator=None, pickled=None):
      if not TRC.init:
         print "WARNING: class Transmission Rows Container not initialized" 
         return None

      self.op_char = op_char
      self.label = label
      self.color = color

      if not operator:
         self.operator = scat.inflate_small_mat(np.identity(TRC.refmodes, dtype="complex"), TRC.modes)
      else:
         self.operator = operator
    
      if not pickled:
         self.eigval, self.eigvec. self.Leigvec = trans.calc_eigensystem(self.operator, TRC.modes)
      elif pickled:
         self.eigvec, self.Leigvec = trans.unpickle_states(pickled, modes, refmodes)

      self.times, self.tisdv = trans.calc_qexpvals(self.eigvec, TRC.q, TRC.qLeigbas,
                                                   TRC.modes, TRC.refmodes, TRC.Nk)
      self.vectransmit, self.vecLtransmit = trans.calc_transmitted_states(TRC.t, TRC.tinv, self.eigvec, self.Leigvec,
                                                                          TRC.modes, TRC.refmodes, TRC.refpos, TRC.Nk)
    
         






qeigval, qeigvec, qLeigvec = trans.calc_eigensystem(qref, modes)
teigval, teigvec, tLeigvec = trans.calc_eigensystem(tref, modes)
modesvec = scat.inflate_small_mat(np.identity(refmodes, dtype="complex"), modes)
qyeigvec, qyLeigvec = trans.unpickle_states(pic_direc+'qpeaks.1.p', modes, refmodes)
#trans.print_freq_indep(qeigval, kmin, kmax, Dk)

eigvec = (qeigvec, teigvec, modesvec, qyeigvec)
Leigvec = (qLeigvec, tLeigvec, modesvec, qyLeigvec)

qeigbas, qLeigbas = trans.calc_qeigbas(q, Nk, modes, refmodes)

Nop = np.shape(eigvec)[0]
times = np.empty((Nop,refmodes,Nk), dtype="complex")
tisdv = np.empty((Nop,refmodes,Nk), dtype="complex")
for i in range(Nop):
    times[i], tisdv[i] = trans.calc_qexpvals(
        eigvec[i], q, qLeigbas, modes, refmodes, Nk)

vectransmit = np.zeros((Nop,Nk,refmodes,modes), dtype="complex")
vecLtransmit = np.zeros((Nop,refmodes,modes), dtype="complex")
for i in range(Nop):
    vectransmit[i], vecLtransmit[i] = trans.calc_transmitted_states(
        t, tinv, eigvec[i], Leigvec[i],
        modes, refmodes, refpos, Nk)


##### output #####

taumax = np.max(np.absolute(qeigval))

labels = (r'$q=-it^{-1}dt/d\omega$',
          r'$t$',
          r'modes',
          r'$\tilde{q}_y$')
color_chars = ('r', 'g', 'b', 'm')
op_chars = ('q','t','m','qy')
op_pretty_chars = (r'$q$',r'$t$',r'$m$',r'$\tilde{q}_y$')


### measure-files (meas) ###

for lbl, opc, oppc, cc, vec, vecL, time, sdv in \
        zip(labels, op_chars, op_pretty_chars, color_chars, vectransmit, vecLtransmit, times, tisdv):
    for state in range(refmodes):
 
        plt.figure(figsize=(16,9))
        plt.suptitle(r'$\rm{Re}(\tau(\omega))$, $\rm{Im}(\tau(\omega))$, and measure $\mu(\omega)$ for '+oppc+'-eigenstates',
                     fontsize=18) 

        plt.subplot2grid((4,100),(0,0), colspan=25)
        height = np.max(np.absolute(time[state])/taumax)
        plt.errorbar( range(Nk), (time[state].real)/taumax, (sdv[state].real)/taumax, fmt=cc+'--d', label=lbl)

        omegaref = plt.vlines(refpos, 0.0, 1.5, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        plt.ylabel(r'$\rm{Re}(\tau(\omega))/|\eta_{\,\rm{max}}(\omega_0)|$')
        plt.ylim(0.0,1.5)


        plt.subplot2grid((4,100),(0,36), colspan=25)
        #height = np.max(np.absolute(time[state].imag)/taumax)
        height = 0.05
        plt.plot( range(Nk), np.absolute(time[state].imag)/taumax, cc+'--d', label=lbl)
        
        plt.vlines(refpos, 0., height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        plt.ylabel(r'$\rm{Im}(\tau(\omega))/|\eta_{\,\rm{max}}(\omega_0)|$')
        plt.ylim(0.,height)


        ax = plt.subplot2grid((4,100),(0,74), colspan=25)
      
        
        mu = np.absolute([1-abs(np.dot(vecL[state], vec[w][state])) for w in range(Nk)])
        #height = np.max(mu)
        height = 0.50
        
        ax.plot(mu, cc+'--d')
        ax.plot([0. for w in range(Nk)])

        ax.vlines(refpos, 0., height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        plt.ylabel(r'$\mu(\omega)$')
        plt.ylim(0.0,height)
        ax.legend(bbox_to_anchor=(1.4,1.1))

    
        plt.subplot2grid((4,100),(2,0), rowspan=2, colspan=100)
        plt.title(r'$\left| \Delta q \left(\omega_0\right) / \eta_{\,\rm{max}}(\omega_0) \right|=$'
                  '%.5f' % np.absolute(sdv[state][refpos]/taumax))
        try:
            plt.tick_params(which='both', bottom='off', top='off', left='off',
                            right='off', labelbottom='off', labelleft='off')
            if (opc == 'q'):
                image = plt.imread(pic_direc+"qeigenstates/"+\
                                       "pic."+pic_filen+".0002.%04i.noise.a1.layer0.jpg" % state)
            elif (opc == 't'):
                image = plt.imread(pic_direc+"teigenstates/"+\
                                       "pic."+pic_filen+".%04i.noise.a1.layer0.jpg" % state)
            elif (opc == 'm'):
                image = plt.imread(pic_direc+"meigenstates/"+\
                                       "pic."+pic_filen+".0000.%04i.streu.layer0.jpg" % state)
            elif (opc == 'qy'):
                image = plt.imread(pic_direc+"peakwavefuncs/"+\
                                       "pic."+pic_filen+".qpeaks.1.%04i.coeff.layer0.jpg" % state)
               
            plt.imshow(image, origin='left')
        except:
            print "WARNING: unable to insert image "+opc+".%i" % state

        plt.savefig(plot_direc+"times."+opc+".%i.jpg" % state)
        plt.clf()


### times-files (times) ###



if (fullout==True):
#####
#
### transmitted coefficients-files (transmit) ##
#
#####
    for state in range(modes):
        height = np.max( np.absolute( np.take( qtransmit, [state], axis=1) ) )
        for i in range(modes):
            plot(\
                np.absolute( np.take( np.take( qtransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
                    '--d' )
            vlines(refpos, 0., 2*height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
            title('absolute values of (lead-mode-)coefficients $c$ at exit')
            xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
            ylabel(r'$|c(\omega)|$')
            legend( bbox_to_anchor=(1.,1.) )
            savefig(direc+"qtransmit.%i.jpg"%state)
            clf()

    for state in range(modes):
        height = np.max( np.absolute( np.take( Qtransmit, [state], axis=1) ) )
        for i in range(modes):
            plot(\
                np.absolute( np.take( np.take( Qtransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
                    '--d' )
            vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
            title('absolute values of (lead-mode-)coefficients $c$ at exit')
            xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
            ylabel(r'$|c(\omega)|$')
            legend( bbox_to_anchor=(1.,1.) )
            savefig(direc+"Qtransmit.%i.jpg"%state)
            clf()
                    
    for state in range(modes):
        height = np.max( np.absolute( np.take( ttransmit, [state], axis=1) ) )
        for i in range(modes):
            plot(\
            np.absolute( np.take( np.take( ttransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
                '--d' )
            vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
            title('absolute values of (lead-mode-)coefficients $c$ at exit')
            xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
            ylabel(r'$|c(\omega)|$')
            legend( bbox_to_anchor=(1.,1.) )
            savefig(direc+"ttransmit.%i.jpg"%state)
            clf()

    for state in range(modes):
        height = np.max( np.absolute( np.take( mtransmit, [state], axis=1) ) )
        for i in range(modes):
            plot(\
                np.absolute( np.take( np.take( mtransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
                    '--d' )
            vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
            title('absolute values of (lead-mode-)coefficients $c$ at exit')
            xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
            ylabel(r'$|c(\omega)|$')
            legend( bbox_to_anchor=(1.,1.) )
            savefig(direc+"mtransmit.%i.jpg"%state)
            clf()


