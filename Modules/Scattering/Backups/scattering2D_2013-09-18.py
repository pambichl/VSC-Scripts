import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from cmath import *

I = 1.J
Pi = np.pi


def read_S(direc, filen, modes_max, returndims=0):
    Sstr  = open(direc+"Smat."+filen+".dat", "r").readlines()
    NS = int(Sstr[0])

    modes = np.zeros((NS))
    mod_indx = -1

    S = np.zeros((NS,2*modes_max,2*modes_max), dtype='complex')
    S_indx = -1

    for line in Sstr[1:-1]:
        line = line.split()

        if len(line) == 1:
            mod_indx += 1
            S_indx   += 1
            if S_indx >= NS: break
            modes[mod_indx] = int(line[0])
        if len(line) == 4:
            row = int(line[0])
            col = int(line[1])
            if row >= modes[mod_indx]/2: row += modes_max - modes[mod_indx]/2 
            if col >= modes[mod_indx]/2: col += modes_max - modes[mod_indx]/2
            val = float(line[2]) + I*float(line[3])
            S[S_indx][row][col] = val
            
    if returndims==0:
        return S
    else:
        return S, modes


def calc_Q(S, dk, inv=False):
    NQ = len(S)/3
    if (len(S)%3 != 0): print "CALC_Q: WARNING, number of S matrices not divisible by 3"

    m = np.shape(S)[1]
    Q = np.zeros((NQ,m,m), dtype='complex')

    for i in range(NQ):
        if (inv==True): 
            Q[i] = -I * np.mat(np.linalg.inv(S[3*i+1])) * np.mat(S[3*i+2] - S[3*i+0]) / (2*dk)
        else:
            Q[i] = -I * np.mat(S[3*i+1]).conj().transpose() * np.mat(S[3*i+2] - S[3*i+0]) / (2*dk)
        Q[i] = np.mat(Q[i])

    return Q


def inflate_vec(v, dim):
    if (len(np.shape(v))!=1):
        print "INFLATE_VEC: BREAK, no vector object passed"
        return
    vinfl = np.zeros((dim), dtype="complex")
    n = np.shape(v)[0]
    vinfl[0:n] = v

    return vinfl


def inflate_small_mat(A, dim):
    if (len(np.shape(A))!=2):
        print "INFLATE_MAT: BREAK, no matrix object passed"
        return
    n,m = np.shape(A)
    if (n!=m):
        print "INFLATE_MAT: BREAK, passed matrix not square matrix"
        return
    Ainfl = np.zeros((dim,dim), dtype="complex")
    Ainfl[0:n,0:m] = A

    return Ainfl


def inflate_large_mat(A, dim):
    if (len(np.shape(A))!=2):
        print "INFLATE_MAT: BREAK, no matrix object passed"
        return
    if (np.shape(A)[0]!=np.shape(A)[1]):
        print "INFLATE_MAT: BREAK, passed matrix not square matrix"
        return 
    if (dim%2!=0):
        print "INFLATE_MAT: BREAK, dimension of inflated matrix not divisible by 2"
        return 
    Ainfl = np.zeros((dim,dim), dtype="complex")
    n = np.shape(A)[0]/2
    m = int(dim/2)

    print n,m
    print range(-m,-m+n+1)
    print A[-n:,-n:]

    Ainfl[0:n,0:n]         = A[0:n,0:n]
    Ainfl[0:n,-m:-m+n]     = A[0:n,-n:]
    Ainfl[-m:-m+n,0:n]     = A[-n:,0:n]
    Ainfl[-m:-m+n,-m:-m+n] = A[-n:,-n:]

    return Ainfl


def assemble_QT(Qlist):
    dim      = np.shape(Qlist)[0]
    blockdim = np.shape(Qlist)[1]

    QT = np.empty([dim*blockdim,dim*blockdim], dtype='complex')
    for i in range(dim):
        for j in range(dim):
            for mode_i in range(blockdim):
                for mode_j in range(blockdim):
                    if (i>=j): QT[i*blockdim + mode_i, j*blockdim + mode_j] = Qlist[i-j][mode_i][mode_j]
                    if (i< j): QT[i*blockdim + mode_i, j*blockdim + mode_j] = (Qlist[j-i].conj().transpose())[mode_i][mode_j]

    return np.mat(QT)


def order_QTeigvec(eigvec, modes):
    dim = np.shape(eigvec)[1]    
    tsteps = int(dim/modes)

    ord_eigvec = np.empty((dim,modes,tsteps), dtype="complex")

    for i in range(dim):
        for j in range(modes):
            ord_eigvec[i][j] = np.take(eigvec[i], range(j,dim,modes))

    return np.array(ord_eigvec)
    

def transform_packcoeff(vec, freqs, times, ny, dy):
    dim = np.shape(freqs)[0]
    modes = np.shape(vec)[0]
    length = np.shape(vec)[1]

    dt = times[1] - times[0]

    vecFT = np.zeros((modes,dim), dtype="complex")
    for n in range(modes):
        for w in range(dim):
            for t in range(length):
                if (freqs[w]*dy < (n+1)*Pi/(ny+1)): break
                vecFT[n][w] += dt *\
                    vec[n][t] *\
                    exp(I*freqs[w]*times[t])

    return vecFT


def calc_initial(vec, freqs, times, nx, ny, dy):
    dim = np.shape(freqs)[0]
    modes = np.shape(vec)[0]
    length = np.shape(vec)[1]
    dx = dy

    dw = freqs[1] - freqs[0]

    chi = np.empty((ny,modes), dtype="complex")
    for i in range(ny):
        for j in range(modes):
            chi[i][j] = sin( (j+1) * Pi / (ny+1) * (i+1) )

    kx = np.empty((modes,dim), dtype="complex")
    for i in range(modes):
        for j in range(dim):
            freqxsqrd = freqs[j]**2 - 2./dy**2 * (1 - cos( (i+1)*Pi/(ny+1) ) )
            kx[i][j] = 1./dx * acos( 1 - freqxsqrd*dx**2/2. )

    v = np.empty((modes,dim), dtype="complex")
    for i in range(modes):
        for j in range(dim):
            v[i][j] = sin(kx[i][j]*dx)/dx

    vecFT = transform_packcoeff(vec, freqs, times, ny, dy)

#    psi_x = np.zeros((modes,2*nx), dtype="complex")
#    for n in range(modes):
#        for ix in range(2*nx):
#            for w in range(dim):
#                for t in range(length):
#                    psi_x[n][ix] +=\
#                        dt*dw/(2.*Pi) *\
#                        1./sqrt(v[n][w]) *\
#                        vec[n][t] *\
#                        exp(I * (kx[n][w]*(ix-nx)*dx + freqs[w]*times[t]) )

    psi_x = np.zeros((modes,2*nx), dtype="complex")

    for n in range(modes):
        for ix in range(2*nx):
            for w in range(dim):
                if (abs(kx[n][w].imag) < 10**-8):
                    psi_x[n][ix] +=\
                        dw/(2.*Pi) *\
                        1./sqrt(v[n][w]) *\
                        vecFT[n][w] *\
                        exp( I*kx[n][w]*(ix-nx)*dx )


    psi = np.zeros((2*nx,ny), dtype="complex")
    for n in range(modes):
        psi += np.dot(np.mat(psi_x[n]).transpose(), np.mat(chi[:,n]))

    return psi
