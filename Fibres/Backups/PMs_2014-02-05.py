#!/usr/bin/env python

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from cmath import *
import matplotlib.pyplot as plt
import sys

from Utils import utils as ut
from Scattering import scattering2D as scat


Pi = np.pi
I = 0. + 1.J


if (len(sys.argv)!=5):
    sys.exit("ABORT: VSC (1,2), parameters filen, pic-filen, # of modes to be considered needed")

machine = int(sys.argv[1]) # on which machine calculation was performed, VSC or VSC2
filen = str(sys.argv[2]) # namestring of calculation where scatter data is stored
pic_filen = str(sys.argv[3]) # namestring of calculation where waveplots are stored
refmodes = int(sys.argv[4]) # number of modes to be considered

if machine == 1:
    data_direc = "../../VSC/Fibres-20130819/" + filen + "/"
elif machine ==2:
    data_direc = "../../VSC2/Fibres-20130802/" + filen + "/"
fullout = False

### needs to be adapted, but output quantities dimensionless! ###
dk = 10**(-8)
### ###

S, dims = scat.read_S(data_direc, filen, old_ver=1)
modes = int(max(dims)/2)

if (len(dims)%3 == 0):
    Nk = len(dims)/3
else:
    sys.exit("ABORT: # of S matrices not divisible by 3")

if (refmodes > int(min(dims)/2)):
    print "WARNING: # of modes under consideration (%i) > # minimum number of modes (%i)" % (refmodes,int(min(dims)/2))


t = []
t_q = []
for i in range(3*Nk):
    m = dims[i]/2
    n = min(m, refmodes)
    t.append( S[i][m:,0:m] )      # full t-matrix
    t_q.append( S[i][m:m+n,0:n] ) # considered part of t-matrix 

q = np.array([scat.calc_Q([t_q[3*i+0],t_q[3*i+1],t_q[3*i+2]], dk, inv=True)[0] for i in range(Nk)])
refpos = int(0.5*(Nk-1)) # k[refpos] = kmean for odd Nk


### define operators at reference frequency ### 
qref = q[refpos]
tref = t_q[3*refpos+1]
tinv = np.linalg.inv(tref)

### calculate and sort eigenvalues and left and right eigenvectors of operators ###
qeigval, qeigvec = ut.sort_eig(qref)
teigval, teigvec = ut.sort_eig(tref)

qeigvec  = qeigvec.transpose()
teigvec  = teigvec.transpose()
modesvec = np.identity(refmodes, dtype="complex")

teigval = ut.sort_eig(t_q[3*(refpos)+1])[0]
np.savetxt(data_direc+"teigvals."+filen+".dat", np.array([range(np.shape(teigval)[0]), teigval.real, teigval.imag]).transpose())

# left eigenvectors are to be calculated before inflating, needed for measure \mu
qLeigvec = np.linalg.inv(qeigvec).transpose()
tLeigvec = np.linalg.inv(teigvec).transpose()

### inflate eigenvectors to maximum size and number by appending zeros ###
qeigvec  = scat.inflate_small_mat(qeigvec, modes)
teigvec  = scat.inflate_small_mat(teigvec, modes)
modesvec = scat.inflate_small_mat(modesvec, modes)

qLeigvec = scat.inflate_small_mat(qLeigvec, modes)
tLeigvec = scat.inflate_small_mat(tLeigvec, modes)


# taumin = np.min( qeigval.real )
# taumax = np.max( qeigval.real )
# Dk = (kmax-kmin)/(Nk-1)
# print 
# print 'estimate for smallest frequency-independence'
# print 'wmin =', kmin, 'wmax =', kmax, 'dw =', Dk
# print 'tau_min =', taumin
# print 'tau_max =', taumax
# print 'frequency range: ', 1./(taumax-taumin)
# print 'in units of dw: ', 1./(taumax-taumin)/Dk
# print


### calculate left and right eigenbasis of q at each frequency, needed for left vectors ###
qeigbas  = []
qLeigbas = []
for i in range(Nk):
    evc = ut.sort_eig(q[i])[1].transpose()
    qeigbas.append(scat.inflate_small_mat(evc, modes)[0:refmodes])
    qLeigbas.append(scat.inflate_small_mat(np.linalg.inv(evc).transpose(), modes)[0:refmodes])


### calculate expectation values of q with proper inner product at each frequency ###
times = np.empty( (3,modes,Nk), dtype="complex")
tisqr = np.empty( (3,modes,Nk), dtype="complex")
tisdv = np.empty( (3,modes,Nk), dtype="complex")
eigvec = (qeigvec, teigvec, modesvec)
for w in range(Nk):
    for state in range(refmodes):
        for op in range(3):
            vec = scat.inflate_vec(eigvec[op][state], modes)
            A   = scat.inflate_small_mat(q[w], modes)

            vecL = 1./np.linalg.norm(np.dot(qLeigbas[w], vec))**2 *\
                np.dot(qLeigbas[w].transpose(), np.dot(qLeigbas[w].conj(), vec.conj()))

            times[op][state][w] = np.dot( vecL, np.dot( A, vec ) )
            tisqr[op][state][w] = np.dot( vecL, np.dot( A, np.dot( A, vec ) ) )
            tisdv[op][state][w] = sqrt( tisqr[op][state][w] - times[op][state][w]**2 )


qtransmit = np.zeros( (Nk,refmodes,modes), dtype="complex" )
ttransmit = np.zeros( (Nk,refmodes,modes), dtype="complex" )
mtransmit = np.zeros( (Nk,refmodes,modes), dtype="complex" )


### calculating output states for unchanged input at different frequencies ###
# taking here the full t-matrices gives at \omega_0 also a measure for scattering into
# unconsidered modes
for i in range(Nk):
    for j in range(refmodes):
        #t_inf = scat.inflate_small_mat(t_q[i], modes)
        t_inf = scat.inflate_small_mat(t[3*i+1], modes)
        qtransmit[i][j] = np.dot(t_inf, qeigvec[j]) / np.linalg.norm(np.dot(t_inf, qeigvec[j]))
        ttransmit[i][j] = np.dot(t_inf, teigvec[j]) / np.linalg.norm(np.dot(t_inf, teigvec[j]))
        mtransmit[i][j] = np.dot(t_inf, modesvec[j]) / np.linalg.norm(np.dot(t_inf, modesvec[j]))

# # taking only large frequency steps (not small ones for derivative)
# qtransmit = np.take(qtransmit, np.arange(1,3*Nk,3), axis=0)
# ttransmit = np.take(ttransmit, np.arange(1,3*Nk,3), axis=0)
# mtransmit = np.take(mtransmit, np.arange(1,3*Nk,3), axis=0)

### left eigenstates times the inverse of the reference t(w_0) ###
qLtransmit = np.zeros( (refmodes,modes), dtype="complex" )
tLtransmit = np.zeros( (refmodes,modes), dtype="complex" )
mLtransmit = np.zeros( (refmodes,modes), dtype="complex" )

tfull_ref_inf = scat.inflate_small_mat(t[3*refpos], modes)
tinv_inf = scat.inflate_small_mat(tinv, modes)
for j in range(refmodes): 
    qLtransmit[j] = np.dot(qLeigvec[j], tinv_inf) * np.linalg.norm(np.dot(tfull_ref_inf, qeigvec[j]))
    tLtransmit[j] = np.dot(tLeigvec[j], tinv_inf) * np.linalg.norm(np.dot(tfull_ref_inf, teigvec[j]))
    mLtransmit[j] = np.dot(modesvec[j], tinv_inf) * np.linalg.norm(np.dot(tfull_ref_inf, modesvec[j]))


print np.dot(qLtransmit[0], qtransmit[refpos,0])



##### output #####

taumax = np.max(np.absolute(qeigval))
plot_direc = data_direc
if machine == 1:
    pic_direc = "../../VSC/Fibres-20130819/" + pic_filen + "/"
elif machine ==2:
    pic_direc = "../../VSC2/Fibres-20130802/" + pic_filen + "/"

labels =\
    (r'$q=-it^{-1}dt/d\omega$',\
    r'$t$',\
    r'modes')
color_chars = ('r', 'g', 'b')
op_chars = ('q','t','m')

vecLtransmit = (qLtransmit, tLtransmit, mLtransmit)
vectransmit  = (qtransmit, ttransmit, mtransmit)


### measure-files (meas) ###

for lbl, opc, cc, vec, vecL, time, sdv in \
        zip(labels, op_chars, color_chars, vectransmit, vecLtransmit, times, tisdv):
    for state in range(refmodes):
 
        plt.figure(figsize=(16,9))
        plt.suptitle(r'$\rm{Re}(\tau(\omega))$, $\rm{Im}(\tau(\omega))$, and measure $\mu(\omega)$ for '+opc+'-eigenstates',
                     fontsize=18) 

        plt.subplot2grid((4,100),(0,0), colspan=25)
        height = np.max(np.absolute(time[state])/taumax)
        plt.errorbar( range(Nk), (time[state].real)/taumax, (sdv[state].real)/taumax, fmt=cc+'--d', label=lbl)

        omegaref = plt.vlines(refpos, 0., 1., colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        plt.ylabel(r'$\rm{Re}(\tau(\omega))/|\eta_{\,\rm{max}}(\omega_0)|$')


        plt.subplot2grid((4,100),(0,36), colspan=25)
        #height = np.max(np.absolute(time[state].imag)/taumax)
        height = 0.05
        plt.plot( range(Nk), np.absolute(time[state].imag)/taumax, cc+'--d', label=lbl)
        
        plt.vlines(refpos, 0., height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        plt.xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        plt.ylabel(r'$\rm{Im}(\tau(\omega))/|\eta_{\,\rm{max}}(\omega_0)|$')
        plt.ylim(0.,0.05)


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
        try:
            plt.tick_params(which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')
            if (opc == 'q'):
                image = plt.imread(pic_direc+opc+"eigenstates/"+"pic."+pic_filen+".0002.%04i.noise.a1.layer0.jpg" % state)
            elif (opc == 't'):
                image = plt.imread(pic_direc+opc+"eigenstates/"+"pic."+pic_filen+".%04i.noise.a1.layer0.jpg" % state)
            elif (opc == 'm'):
                image = plt.imread(pic_direc+opc+"eigenstates/"+"pic."+pic_filen+".0000.%04i.streu.layer0.jpg" % state)

                plt.imshow(image, origin='left')
        except Exception:
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


