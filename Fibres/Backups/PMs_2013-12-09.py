#!/usr/bin/env python

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from cmath import *
from matplotlib.pyplot import *
import sys

from Utils import utils as ut
from Scattering import scattering2D as scat


print "SHJDSHJDHJDHJDH"

Pi = np.pi
I = 0. + 1.J


if (len(sys.argv)!=3):
    sys.exit("ABORT: parameters filen, # of radius variations, # of configurations,  # of frequency variations needed")

filen = str(sys.argv[1]) # namestring of calculation
refmodes = int(sys.argv[2]) # number of modes to be considered

data_direc = "../../VSC2/Fibres-20130802/" + filen + "/scatterdata/"
fullout = False

#modes_min = 3.40
#modes_max = 3.60

#modes=int(modes_max)


#ppwhl = 10

#W = 1.
#L = 10.*W

#kmin  = modes_min * Pi / W
#kmax  = modes_max * Pi / W
#kmean = 0.5*(kmax+kmin)
#dk = kmean * 10**-8
#Nk = 7

dk = 0.00023626

S, dims = scat.read_S(direc, filen)
modes = int(max(dims)/2)

if (refmodes > int(min(dims)/2):
        sys.exit("ABORT: # of modes under consideration ($i) > # number of modes (%i)" % (refmodes,modes))

exit

t = []
t_q = []
for i in range(3*Nk):
    m = dimlist[i]/2
    n = min(m, refmodes)
    t.append( S[i][m:,0:m] )      # full t-matrix
    t_q.append( S[i][m:m+n,0:n] ) # considered part of t-matrix 

q = []
for i in range(Nk):
    q.append( scat.calc_Q([t_q[3*i+0],t_q[3*i+1],t_q[3*i+2]], dk, inv=True)[0] )


refpos = int(0.5*(Nk-1)) # k[refpos] = kmean for odd Nk


### define operators at reference frequency ### 
qref = q[refpos]
tref = t_q[3*refpos+1]
tinv = np.linalg.inv(tref)

### calculate and sort eigenvalues and left and right eigenvectors of operators ###
qeigval = ut.sort_eig(qref)[0]
teigval = ut.sort_eig(tref)[0]
#teigval = np.linalg.eig(tref)[0]

teigval = ut.sort_eig(t_q[3*(refpos)+0])[0]
np.savetxt(direc+"teigvals."+filen+".dat", np.array([range(np.shape(teigval)[0]), teigval.real, teigval.imag]).transpose())

qeigvec  = ut.sort_eig(qref)[1].transpose()
teigvec  = ut.sort_eig(tref)[1].transpose()
#teigvec  = np.linalg.eig(tref)[1].transpose()
modesvec = np.identity(refmodes, dtype="complex")

# left eigenvectors are to be calculated before inflating
qLeigvec = np.linalg.inv(qeigvec).transpose()
tLeigvec = np.linalg.inv(teigvec).transpose()


### inflate eigenvectors to maximum size and number by appending zeros ###
qeigvec  = scat.inflate_small_mat(qeigvec, modes)
teigvec  = scat.inflate_small_mat(teigvec, modes)
modesvec = scat.inflate_small_mat(modesvec, modes)

qLeigvec = scat.inflate_small_mat(qLeigvec, modes)
tLeigvec = scat.inflate_small_mat(tLeigvec, modes)

eigvec = (qeigvec, teigvec, modesvec)


taumin = np.min( qeigval.real )
taumax = np.max( qeigval.real )
Dk = (kmax-kmin)/(Nk-1)
print 
print 'estimate for smallest frequency-independence'
print 'wmin =', kmin, 'wmax =', kmax, 'dw =', Dk
print 'tau_min =', taumin
print 'tau_max =', taumax
print 'frequency range: ', 1./(taumax-taumin)
print 'in units of dw: ', 1./(taumax-taumin)/Dk
print


#Qt = scat.calc_Q(t,dk,inv=False)
#print np.linalg.eig(Qt[refpos])[0]
#for w in range(Nk): Qt[w] = 0.5 * ( Qt[w] + Qt[w].conj().transpose() )
#print np.linalg.eig(Qt[refpos])[0]
#print np.linalg.eig(Q[refpos])[0]


qeigbas  = []
qLeigbas = []
for i in range(Nk):
    evc = ut.sort_eig(q[i])[1].transpose()
    qeigbas.append(scat.inflate_small_mat(evc, modes)[0:refmodes])
    qLeigbas.append(scat.inflate_small_mat(np.linalg.inv(evc).transpose(), modes)[0:refmodes])


times = np.empty( (3,modes,Nk), dtype="complex")
tisqr = np.empty( (3,modes,Nk), dtype="complex")
tisdv = np.empty( (3,modes,Nk), dtype="complex")
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


# <codecell>

qtransmit = np.zeros( (3*Nk,refmodes,modes), dtype="complex" )
Qtransmit = np.zeros( (3*Nk,refmodes,modes), dtype="complex" )
ttransmit = np.zeros( (3*Nk,refmodes,modes), dtype="complex" )
mtransmit = np.zeros( (3*Nk,refmodes,modes), dtype="complex" )

### calculating output states for unchanged input at different frequencies ###
# taking here the full t-matrices gives at w_0 a measure for scattering into
# unconsidered modes
for i in range(3*Nk):
    for j in range(refmodes):
        #t_inf = scat.inflate_small_mat(t_q[i], modes)
        t_inf = scat.inflate_small_mat(t[i], modes)
        qtransmit[i][j] = np.dot(t_inf, qeigvec[j]) / np.linalg.norm(np.dot(t_inf, qeigvec[j]))
        ttransmit[i][j] = np.dot(t_inf, teigvec[j]) / np.linalg.norm(np.dot(t_inf, teigvec[j]))
        mtransmit[i][j] = np.dot(t_inf, modesvec[j]) / np.linalg.norm(np.dot(t_inf, modesvec[j]))

# taking only large frequency steps (not small ones for derivative)
qtransmit = np.take(qtransmit, np.arange(1,3*Nk,3), axis=0)
ttransmit = np.take(ttransmit, np.arange(1,3*Nk,3), axis=0)
mtransmit = np.take(mtransmit, np.arange(1,3*Nk,3), axis=0)


### left eigenstates times the inverse of the reference t(w_0) ###
qLtransmit = np.zeros( (refmodes,modes), dtype="complex" )
tLtransmit = np.zeros( (refmodes,modes), dtype="complex" )
mLtransmit = np.zeros( (refmodes,modes), dtype="complex" )

tref_inf = scat.inflate_small_mat(tref, modes)
tinv_inf = scat.inflate_small_mat(tinv, modes)
for j in range(refmodes): 
    qLtransmit[j] = np.dot(qLeigvec[j], tinv_inf) * np.linalg.norm(np.dot(tref_inf, qeigvec[j]))
    tLtransmit[j] = np.dot(tLeigvec[j], tinv_inf) * np.linalg.norm(np.dot(tref_inf, teigvec[j]))
    mLtransmit[j] = np.dot(modesvec[j], tinv_inf) * np.linalg.norm(np.dot(tref_inf, modesvec[j]))

# <codecell>

print np.dot(qLtransmit[0], qtransmit[refpos,0])
#print qtransmit[0][1]
#print np.dot( tinv, tref )
#print np.dot( tinv_inf, tref_inf )

#(np.dot(t_inf, qeigvec[j]



##### output #####



labels =\
    (r'$q=-it^{-1}dt/d\omega$',\
    r'$t$',\
    r'modes')
colorchar = ('r', 'g', 'b')


### measure-files (meas) ###

for op in range(3):
    for state in range(refmodes):
        if (op==0):
            height = np.max(np.absolute(np.array([1-abs(np.dot(qLtransmit[state], qtransmit[i][state])) for i in range(Nk)])))
            plot(\
                np.absolute(np.array(\
                        [1-abs(np.dot(qLtransmit[state], qtransmit[i][state])) for i in range(Nk)]\
                            )),\
                    colorchar[op]+'--d', label=labels[op])
        elif (op==1):
            height = np.max(np.absolute(np.array([1-abs(np.dot(tLtransmit[state], ttransmit[i][state])) for i in range(Nk)])))
            plot(\
                np.absolute(np.array(\
                        [1-abs(np.dot(tLtransmit[state], ttransmit[i][state])) for i in range(Nk)]\
                            )),\
                    colorchar[op]+'--d', label=labels[op])
        elif (op==2):
            height = np.max(np.absolute(np.array([1-abs(np.dot(mLtransmit[state], mtransmit[i][state])) for i in range(Nk)])))
            plot(\
                np.absolute(np.array(\
                        [1-abs(np.dot(mLtransmit[state].transpose(), mtransmit[i][state])) for i in range(Nk)]\
                            )),\
                    colorchar[op]+'--d', label=labels[op])

        plot([0. for i in range(Nk)])
        vlines(refpos, 0., height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        title(r'measure $\mu$ for $\omega$-dependance of output states')
        xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        ylabel(r'$\mu(\omega)$')
        ylim(0.0,0.10)
        legend( bbox_to_anchor=(1.,1.) )
        savefig(direc+"meas.%i.%i.jpg"%(op,state))
        clf()


### times-files (times) ###


for op in range(3):
    for state in range(refmodes):
        height = np.max( np.absolute(times[op][state]) )
        #plot(times[op][state], 'r--d', label=labels[op])
        errorbar( range(Nk), times[op][state].real, tisdv[op][state].real, fmt=colorchar[op]+'--d', label=labels[op])

        vlines(refpos, 0.75*height, 1.25*height, colors='yellow', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
        title(r'$q(\omega_0)$-expectation values for '+labels[op]+'-eigenstates')
        xlabel(r'$(\omega-\omega_{min})/\Delta\omega$')
        ylabel(r'$\tau(\omega)$')
        #ylim(0.5*height,1.5*height)
        ylim(taumin,taumax)
        legend( bbox_to_anchor=(1.,1.) )
        savefig(direc+"times.%i.%i.jpg"%(op,state))
        clf()


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
            vlines(refpos, 0., height, colors='red', linestyles='--', linewidth=2.0, label=r'$\omega_0$')
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


