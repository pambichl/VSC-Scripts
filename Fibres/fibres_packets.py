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

opcList = ('q', 't', 'm', 'qy', 'tdty', 'Kfull', 'Kcent')
if (len(sys.argv) > 1):
    # which states are to be plotted
   opc = str(sys.argv[1])
if ( (len(sys.argv) <= 1) or (opc not in opcList) ): 
   print "choose states to be plotted:"
   print "q... q-eigenstates"
   print "t... t-eigenstates"
   print "m... modes"
   print "qy... q_y-eigenstates"
   print "tdty... t_y^dagger.t_y-eigenstates"
   print "Kfull... K-eigenstates for full frequency range"
   print "Kcent... K-eigenstates for center frequencies"
   print "WARNING: assuming modes to be desired states"
   print
   opc = 'm'

fullout = True

try:
    # on which machine calculation was performed, VSC or VSC2
    machine = int(sys.argv[2])
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


class TRC:
   '''
   Transmission Rows Container, class that calculates and
   subsequently contains all data needed for output.
   Needs to be intitialized.
   If no operator is given and no pickled vectors are
   read in, modes are assumed to be the desired states.
   '''
   init = False

   @classmethod
   def class_init(cls, modes, refmodes, refpos,
                  Nk, t, tref, q):

      TRC.modes = modes
      TRC.refmodes = refmodes
      TRC.refpos = refpos
      TRC.Nk = Nk
      TRC.t = t
      TRC.tref = tref
      TRC.tinv = np.linalg.inv(tref)
      TRC.q = q
      TRC.init = True

      TRC.qeigbas, TRC.qLeigbas = trans.calc_qeigbas(q, Nk, modes, refmodes)

   def __init__(self, op_char='m', label=r'modes', color='b',
                operator=np.array(None), pickled=None):

      if not TRC.init:
         print "WARNING: class Transmission Rows Container not initialized" 
         return None

      self.op_char = op_char
      self.label = label
      self.color = color

      if not operator.any() and not pickled:
         self.operator = np.identity(TRC.refmodes, dtype="complex")
         self.eigvec = scat.inflate_small_mat(np.identity(TRC.refmodes, dtype="complex"), TRC.modes)
         self.Leigvec = scat.inflate_small_mat(np.identity(TRC.refmodes, dtype="complex"), TRC.modes)
      elif not operator.any() and pickled:
         self.operator = np.identity(TRC.refmodes, dtype="complex")
         self.eigvec, self.Leigvec = trans.unpickle_states(pickled, modes, refmodes)
      else:
         self.operator = operator
         self.eigval, self.eigvec, self.Leigvec = trans.calc_eigensystem(self.operator, TRC.modes)

      self.times, self.tisdv = trans.calc_qexpvals(self.eigvec, TRC.q, TRC.qLeigbas,
                                                   TRC.modes, TRC.refmodes, TRC.Nk)
      self.vectransmit, self.vecLtransmit = trans.calc_transmitted_states(TRC.t, TRC.tinv, self.eigvec, self.Leigvec,
                                                                          TRC.modes, TRC.refmodes, TRC.refpos, TRC.Nk)

      self.mu = trans.calc_mu(self.vectransmit, self.vecLtransmit, TRC.refmodes, TRC.Nk)

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

taumax = np.max(np.absolute(np.linalg.eig(qref)[0]))
        
TRC.class_init(modes, refmodes, refpos, Nk, t, tref, q)
qCont = TRC('q', r'$q$', 'r', operator=qref)
tCont = TRC('t', r'$t$', 'g', operator=tref)
mCont = TRC('m', 'modes', 'b')
qyCont = TRC('qy', r'$\tilde{q}_y$', 'm', pickled=pic_direc+'qpeaks.1.p')
tdtyCont = TRC('tdty', r'$\tilde{t}_y^{\dagger} \tilde{t}_y$', 'c', pickled=pic_direc+'tdtpeaks.1.0.p')
KfCont = TRC('Kfull', r'$K$', 'orange', pickled=data_direc+'KEigVecs.0.p')
KcCont = TRC('Kcent', r'$K$', 'orange', pickled=data_direc+'KEigVecs.2.p')

Cont = (qCont, tCont, mCont, qyCont, tdtyCont, KfCont, KcCont)
ContDict = dict(zip(opcList, Cont))


##### output #####


cont = ContDict[opc]
lbl = cont.label
#opc = cont.op_char
cc = cont.color
vec = cont.vectransmit
vecL = cont.vecLtransmit
time = cont.times
sdv = cont.tisdv
mu = cont.mu
   
### overview plots ###
for state in range(refmodes):
 
   plt.figure(figsize=(16,9))
   plt.suptitle(r'$\rm{Re}(\tau(\omega))$, $\rm{Im}(\tau(\omega))$, and measure $\mu(\omega)$ for '+lbl+'-eigenstates',
                fontsize=18) 

   plt.subplot2grid((4,100),(0,0), colspan=25)
   height = np.max(np.absolute(time[state])/taumax)
   plt.errorbar( range(Nk), (time[state].real)/taumax, (sdv[state].real)/taumax, color=cc, fmt='--d')

   omegaref = plt.vlines(refpos, 0.0, 1.5, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
   plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
   plt.ylabel(r'$\rm{Re}(\tau(\omega))/|\eta_{\,\rm{max}}(\omega_0)|$')
   plt.ylim(0.0,1.5)


   plt.subplot2grid((4,100),(0,36), colspan=25)
        #height = np.max(np.absolute(time[state].imag)/taumax)
   height = 0.05
   plt.plot( range(Nk), np.absolute(time[state].imag)/taumax, '--d', color=cc)
        
   plt.vlines(refpos, 0., height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
   plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
   plt.ylabel(r'$\rm{Im}(\tau(\omega))/|\eta_{\,\rm{max}}(\omega_0)|$')
   plt.ylim(0.,height)


   ax = plt.subplot2grid((4,100),(0,74), colspan=25)
      
        #height = np.max(mu)
   height = 0.50
        
   ax.plot(mu[state], '--d', color=cc)
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
                               "pic."+pic_filen+".0002.%04i.noise.a1.%04i.layer0.jpg" % (state, state+1))
      elif (opc == 't'):
         image = plt.imread(pic_direc+"teigenstates/"+\
                               "pic."+pic_filen+".0002.%04i.noise.a1.%04i.layer0.jpg" % (state, state+1))
      elif (opc == 'm'):
         image = plt.imread(pic_direc+"meigenstates/"+\
                               "pic."+pic_filen+".0000.%04i.streu.%04i.layer0.jpg" % (state, state+1))
      elif (opc == 'qy'):
         image = plt.imread(pic_direc+"peakwavefuncs/"+\
                               "pic."+pic_filen+".qpeaks.1.%04i.coeff.layer0.jpg" % state)
      elif (opc == 'tdty'):
         image = plt.imread(pic_direc+"peakwavefuncs/"+\
                               "pic."+pic_filen+".tdtpeaks.1.0.%04i.coeff.layer0.jpg" % state)
      elif (opc == 'Kfull'):
         image = plt.imread(data_direc+"Keigenstates/"+\
                               "pic."+filen+".full.%04i.coeff.%04i.layer0.jpg" % (state,80*(state+1)))  
      elif (opc == 'Kcent'):
         image = plt.imread(data_direc+"Keigenstates/"+\
                               "pic."+filen+".cent.%04i.coeff.%04i.layer0.jpg" % (state,80*(state+1)))      

      plt.imshow(image, origin='left')
   except:
      print "WARNING: unable to insert image "+opc+".%i" % state

   plt.savefig(plot_direc+"times."+opc+".%i.jpg" % state)
   plt.clf()


### transmitted coefficients plots ###
if (fullout==True):

   renorm_vec = vec
   for w in range(Nk):
      for state in range(refmodes):
         renorm_vec[w,state] = vec[w,state] / np.linalg.norm(vec[w,state,:refmodes].flatten())
   print np.shape(renorm_vec)
   print np.shape(vec)

   for state in range(refmodes):
      height = np.max( np.absolute( np.take( vec, [state], axis=1) ) )
   
      for i in range(modes):
            #plt.plot(np.absolute( np.take( np.take( vec, [state], axis=1), [i], axis=2 ) ).flatten(), '--d', label=r'mode %i'%i)
         plt.plot(np.absolute(renorm_vec[:,state,i]), '--d', label=r'mode %i'%i)
      plt.vlines(refpos, 0., 2*height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
      plt.title('absolute values of (lead-mode-)coefficients $c$ at exit')
      plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
      plt.ylabel(r'$|c(\omega)|$')
      plt.legend( bbox_to_anchor=(1.,1.) )
      plt.savefig(plot_direc+opc+"transmit.%i.jpg"%state)
      plt.clf()

    # for state in range(modes):
    #     height = np.max( np.absolute( np.take( qtransmit, [state], axis=1) ) )
    #     for i in range(modes):
    #         plot(\
    #             np.absolute( np.take( np.take( qtransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
    #                 '--d' )
    #         vlines(refpos, 0., 2*height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
    #         title('absolute values of (lead-mode-)coefficients $c$ at exit')
    #         xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    #         ylabel(r'$|c(\omega)|$')
    #         legend( bbox_to_anchor=(1.,1.) )
    #         savefig(direc+"qtransmit.%i.jpg"%state)
    #         clf()

    # for state in range(modes):
    #     height = np.max( np.absolute( np.take( Qtransmit, [state], axis=1) ) )
    #     for i in range(modes):
    #         plot(\
    #             np.absolute( np.take( np.take( Qtransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
    #                 '--d' )
    #         vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
    #         title('absolute values of (lead-mode-)coefficients $c$ at exit')
    #         xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    #         ylabel(r'$|c(\omega)|$')
    #         legend( bbox_to_anchor=(1.,1.) )
    #         savefig(direc+"Qtransmit.%i.jpg"%state)
    #         clf()
                    
    # for state in range(modes):
    #     height = np.max( np.absolute( np.take( ttransmit, [state], axis=1) ) )
    #     for i in range(modes):
    #         plot(\
    #         np.absolute( np.take( np.take( ttransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
    #             '--d' )
    #         vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
    #         title('absolute values of (lead-mode-)coefficients $c$ at exit')
    #         xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    #         ylabel(r'$|c(\omega)|$')
    #         legend( bbox_to_anchor=(1.,1.) )
    #         savefig(direc+"ttransmit.%i.jpg"%state)
    #         clf()

    # for state in range(modes):
    #     height = np.max( np.absolute( np.take( mtransmit, [state], axis=1) ) )
    #     for i in range(modes):
    #         plot(\
    #             np.absolute( np.take( np.take( mtransmit, [state], axis=1), [i], axis=2 ) ).flatten(),\
    #                 '--d' )
    #         vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
    #         title('absolute values of (lead-mode-)coefficients $c$ at exit')
    #         xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
    #         ylabel(r'$|c(\omega)|$')
    #         legend( bbox_to_anchor=(1.,1.) )
    #         savefig(direc+"mtransmit.%i.jpg"%state)
    #         clf()


