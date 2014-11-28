### functionals ###

# matrices for transverse spread
yvals = np.linspace(0.0,1.0,100)
ky = (np.arange(self.modes)+1)*np.pi
kx = np.array([np.sqrt(w**2 - ky**2) for w in wlist[1::3]])
chi = np.sqrt(2.) * np.einsum('wm,ym->wym', 1.0/np.sqrt(kx),
                              np.sin(np.einsum('m,y->ym', ky, yvals)))
sw_in = np.einsum('w,wym->wym', phi, chi)
sw_out = np.einsum('w,wmn,wyn->wym', phi, t, chi)
self.Y1_in = np.einsum('wym,y,wyn->mn', sw_in.conj(), yvals, sw_in)
self.Y2_in = np.einsum('wym,y,wyn->mn', sw_in.conj(), yvals**2, sw_in)
self.Y1_out = np.einsum('wym,y,wyn->mn', sw_out.conj(), yvals, sw_out)
self.Y2_out = np.einsum('wym,y,wyn->mn', sw_out.conj(), yvals**2, sw_out)

# matrices for temporal correlation
self.C2 = np.einsum('wmn,wml->wnl', dphit.conj(), phit) +\
    np.einsum('wmn,wml->wnl', phit.conj(), dphit)
self.C0 = np.einsum('wmn,wml->wnl', phit.conj(), phit)

# self.wi = [np.einsum('wlm,wln->wmn', self.phit.conj(), self.phit)]
# for o in (1+np.arange(self.wN-1)): 
#     self.wi.append(np.einsum('wlm,wln->wmn', self.phit[:-o].conj(), self.phit[o:]))
# self.wi = np.array(self.wi)

    def G_trans(self, phat):
        '''
        Apply non-linear gradient function for transverse spread
        to vector phat.
        '''

        y1_in = np.dot( phat.conj(), np.dot( self.Y1_in, phat))
        y2_in = np.dot( phat.conj(), np.dot( self.Y2_in, phat))
        y1_out = np.dot( phat.conj(), np.dot( self.Y1_out, phat))
        y2_out = np.dot( phat.conj(), np.dot( self.Y2_out, phat))
        d = np.dot( phat.conj(), np.dot( self.T0, phat))

        G = (-y2_out/d**2 + 2.0*y1_out**2/d**3) * np.dot(self.T0, phat) +\
            1.0/d * np.dot(self.Y2_out, phat) -\
            2.0*y1_out/d**2 * np.dot(self.Y1_out, phat) -\
            ( np.dot(self.Y2_in, phat) -\
                  2.0*y1_in * np.dot(self.Y1_in, phat) )
        
        return G


    def L_trans(self, phat):
        '''
        Functional to be minimized giving transverse spread.
        '''

        y1_in = np.dot( phat.conj(), np.dot( self.Y1_in, phat))
        y2_in = np.dot( phat.conj(), np.dot( self.Y2_in, phat))
        y1_out = np.dot( phat.conj(), np.dot( self.Y1_out, phat))
        y2_out = np.dot( phat.conj(), np.dot( self.Y2_out, phat))
        d = np.dot( phat.conj(), np.dot( self.T0, phat))
        
        L = y2_out/d - (y1_out/d)**2 - ( y2_in - (y1_in)**2 )
        
        return L


    # O0 = np.einsum('w,w', olist**0, R) * self.Dw
    # O1 = np.einsum('w,w', olist**1, R) * self.Dw
    # O2 = np.einsum('w,w', olist**2, R) * self.Dw
    
        ### das Drum sollte ja symmetrisch um Ursprung sein, also EW = 0

        # G = (O2/O0**2 - 2.0*O1**2/O0**3) * np.einsum('w,wn', olist**0, w_mat) * self.Dw -\
        #     1.0/O0 * np.einsum('w,wn', olist**2, w_mat) * self.Dw +\
        #     2.0*O1/O0**2 * np.einsum('w,wn', olist**1, w_mat) * self.Dw

        #Varianz mit eigentlichem EW = 0
        # G = O2/O0**2 * np.einsum('w,wn', olist**0, w_mat) * self.Dw -\
        #     1.0/O0 * np.einsum('w,wn', olist**2, w_mat) * self.Dw

        #EW mit nur positiven o's 
        # G = O1/O0**2 * np.einsum('w,wn', olist**0, w_mat) * self.Dw -\
        #     1.0/O0 * np.einsum('w,wn', olist**1, w_mat) * self.Dw

        # #EW mit nur positiven o's ohne Normierung
        # G = -np.einsum('w,wn', olist**1, w_mat) * self.Dw

        # #Varianz ohne Normierung
        # G = -np.einsum('w,wn', olist**2, w_mat) * self.Dw +\
        #     2.0 * O1 * np.einsum('w,wn', olist**1, w_mat) * self.Dw

        # G = -np.einsum('w,wn', olist**0, w_mat) * self.Dw

        #dR = 1.0/d * np.einsum('w,wn->wn', np.abs(rn+0.0000000001)**(-1.0)/(2.0), w_mat) -\
        #    1.0/d**2 * np.einsum('w,n->wn', np.abs(rn), np.dot(self.T0, phat))
        #dR = 1.0/d * w_mat -\
        #    1.0/d**2 * np.einsum('w,n->wn', np.abs(rn)**2, np.dot(self.T0, phat))
        #dR = w_mat

        #G = np.einsum('w,wn', (olist[1:])**(-1.0/1.0), dR[1:])
        #G = -np.einsum('w,wn', olist, dR)
        #G = -np.einsum('w,wn', olist**2, dR)
        #G = -np.einsum('w,wn', olist**2, dR) - 2.0*np.einsum('w,w,wn', olist, R, dR)
        #G = -np.sum(dR, axis=0)

  # wN = self.wN
        # modes = self.modes
        # wi = self.wi
        # dr = np.zeros((wN,modes), dtype='complex')
        # r = np.zeros((wN), dtype='complex')
        
        # psi = np.einsum('wmn,n->wm', self.phit, phat)

        # ri = np.einsum('wn,wn->w', psi.conj(), psi)
        # r[0] = np.einsum('w,w', ri.conj(), ri) * self.Dw  
        # dr[0] = np.einsum('w,wmn,n->m', ri.conj(), wi[0], phat) * self.Dw +\
        #         np.einsum('w,wnm,n->m', ri, wi[0].conj(), phat) * self.Dw

        # for o in (1+np.arange(wN-1)):
        #     ri = np.einsum('wn,wn->w', psi[:-o].conj(), psi[o:])
        #     r[o] = np.einsum('w,w', ri.conj(), ri) * self.Dw
        #     dr[o] = np.einsum('w,wmn,n->m', ri.conj(), wi[o], phat) * self.Dw +\
        #         np.einsum('w,wnm,n->m', ri, wi[o].conj(), phat) * self.Dw
                     
        # G = -np.sum(dr, axis=0) * self.Dw / r[0] + np.sum(r, axis=0) * self.Dw / r[0]**2 * dr[0]


        #R = np.abs(rn)**2
        #R = np.abs(rn) / d
        #R = np.abs(rn)**2 / d
    
        # On = np.einsum('w,w', (olist[1:])**(-1.0/1.0), R[1:]) * self.Dw
        # O0 = np.einsum('w,w', olist**0, R) * self.Dw
        # O1 = np.einsum('w,w', olist**1, R) * self.Dw
        # O2 = np.einsum('w,w', olist**2, R) * self.Dw

        #L = -O2/O0 + O1**2/O0**2
        #L = -O2/O0
        #L = -O1/O0
        #L = -O1
        #L = -O2
        #L = -O2 - O1**2
        #L = -O0
        #L = O0
        #L = -np.sum(R, axis=0)
        #L = On
        
        # wN = np.shape(self.phit)[0]
        # r = np.zeros((wN), dtype='complex')

        # psi = np.einsum('wmn,n->wm', self.phit, phat)

        # ri = np.einsum('wn,wn->w', psi.conj(), psi)
        # r[0] = np.einsum('w,w', ri.conj(), ri) * self.Dw
           
        # for o in (1+np.arange(wN-1)):
        #     ri = np.einsum('wn,wn->w', psi[:-o].conj(), psi[o:])
        #     r[o] = np.einsum('w,w', ri.conj(), ri) * self.Dw
                     
        # L = -np.sum(r, axis=0) * self.Dw / r[0]

    def G_corr_t(self, phat):
        '''
        Apply non-linear gradient function for spatial stability
        using temporal correlation length to vector phat.
        '''

        h2 = np.einsum('m,wmn,n->w', phat.conj(), self.C2, phat)
        h0 = np.einsum('m,wmn,n->w', phat.conj(), self.C0, phat)
        I2 = np.einsum('w,w', h2, h2)
        I0 = np.einsum('w,w', h0, h0)
        i2 = 2.0 * np.einsum('w,wmn,n->m', h2, self.C2, phat)
        i0 = 2.0 * np.einsum('w,wmn,n->m', h0, self.C0, phat)

        G = I2 / I0**2 * i0 - 1.0 / I0 * i2

        return G 


    def L_corr_t(self, phat):
        '''
        Functional to be minimized for spatial stability
        giving negative square of temporal correlation
        length of output.
        '''

        I2 = np.einsum('m,wmn,n->w', phat.conj(), self.C2, phat) * self.Dw
        I0 = np.einsum('m,wmn,n->w', phat.conj(), self.C0, phat) * self.Dw
        L =  -np.sum(I2**2, axis=0) / np.sum(I0**2, axis=0)

        return L




### opt_movie ###

# temporal autocorrelation function of output signal
def calc_Rc(tau):
   '''
   Calculates absolute square of temporal autocorrelation
   function of calculated output signal. For convenience
   Rc is renormalized by 2 pi.
   '''

   r = np.einsum('wm,wm,w', psi.conj(), psi, np.exp(-I * wlist[1::3] * tau)) * Dw / Nphi
   Rc = np.absolute(r)**2 / (2.0*np.pi)

   r = np.zeros((tN), dtype='complex')
   r[0] = np.einsum('wm,wm', psi_t.conj(), psi_t) * dt
   for tau in (1+np.arange(tN-1)):
      r[tau] =np.einsum('wm,wm', psi_t[:-tau].conj(), psi_t[tau:]) * dt

   return Rc

def calc_rct():
   '''
   Calculates absolute square of temporal autocorrelation
   function of calculated output signal. For convenience
   Rc is renormalized by 2 pi.
   '''

   r = np.zeros((tN), dtype='complex')
   r[0] = np.einsum('wm,wm', psi_t.conj(), psi_t) * dt
   for tau in (1+np.arange(tN-1)):
      r[tau] =np.einsum('wm,wm', psi_t[:-tau].conj(), psi_t[tau:]) * dt

   return r


tau_list = np.arange(tN) * dt
Rct = np.abs(calc_rct())**2 / (2.0*np.pi)
tau0 = np.einsum('t,t', tau_list**0, Rct)
tau1 = np.einsum('t,t', tau_list**1, Rct)
tau2 = np.einsum('t,t', tau_list**2, Rct)

def calc_rcw():
   '''
   Calculates absolute square of temporal autocorrelation
   function of calculated output signal. For convenience
   Rc is renormalized by 2 pi.
   '''

   wN = np.shape(psi)[0]
   r = np.zeros((wN), dtype='complex')
   r[0] = np.einsum('wm,wm', psi.conj(), psi) * Dw
   for o in (1+np.arange(wN-1)):
      r[o] =np.einsum('wm,wm', psi[:-o].conj(), psi[o:]) * Dw

   return r

Ro = np.abs(calc_rcw())**2

olist = np.arange(np.shape(Ro)[0]) * Dw
mean_o = (np.einsum('w,w', olist**1, Ro) * Dw)
#mean_o = (np.einsum('w,w', olist**1, Ro) * Dw) / (np.einsum('w,w', olist**0, Ro) * Dw)

# 1st and 2nd moment of absolute square of temporal autocorrelation function
Rc_t = np.vectorize(calc_Rc)(tlist)
tau = np.einsum('t,t', tlist, Rc_t)
tau_2 = np.einsum('t,t', tlist**2, Rc_t)
vart_corr_out = tau_2 - tau**2

# 1st and 2nd momenta for transverse coordinates of wave packets in time domain
y_in = np.einsum('y,ty,ty', yvals, sigin.conj(), sigin)  / np.sum(np.absolute(sigin)**2)
y2_in = np.einsum('y,ty,ty', yvals**2, sigin.conj(), sigin)  / np.sum(np.absolute(sigin)**2)
y_out = np.einsum('y,ty,ty', yvals, sigout.conj(), sigout)  / np.sum(np.absolute(sigout)**2)
y2_out = np.einsum('y,ty,ty', yvals**2, sigout.conj(), sigout)  / np.sum(np.absolute(sigout)**2)

# variances of transverse coordinates in time domain
vary_in = y2_in - y_in**2
vary_out = y2_out - y_out**2

r"$\left\langle y_{in} \right\rangle = %f$, " % y_in.real +
r"$\sigma_{\perp,in} = %f$, " % np.sqrt(vary_in.real) + 

r"$\left\langle y_{out} \right\rangle = %f$, " % y_out.real +
r"$\sigma_{\perp,out} = %f$, " % np.sqrt(vary_out.real) +
r"$\sigma_{r} = %f$" % np.sqrt(vart_corr_out)

print 
print np.sum(Rc_w/Rc_w[0])
print np.einsum('w,w', olist[1:]**(-1.0), Rc_w[1:])
print np.einsum('w,w', olist, Rc_w)
print np.einsum('w,w', olist**1, Rc_w)
print np.einsum('w,w', olist**2, Rc_w)
print
print np.einsum('w,w', olist[1:]**(-1.0), Rc_w[1:]/Rc_w[0])
print np.einsum('w,w', olist, Rc_w/Rc_w[0])
print np.einsum('w,w', olist**1, Rc_w/Rc_w[0])
print np.einsum('w,w', olist**2, Rc_w/Rc_w[0])
