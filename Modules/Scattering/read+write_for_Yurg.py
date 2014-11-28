import numpy as np

I = 1.J
I = np.complex(0.,1.)
Pi = np.pi


def read_S(direc, filen, old_ver=0):
    '''
    reads in scattering matrices from 'Smat.*.dat' files
    also stores list of leading dimensions of matrices
    and list of scattering energies (=k_0**2/2)

    new version: scattering energies are read in
    old version: scattering energies are not contained in file
    '''
    S_str  = open(direc+"Smat."+filen+".dat", "r").readlines()
    
    data = [map(float, line.strip().split()) for line in S_str if line != '\n']

    S = []
    dim_list = []
    EF_list = []

    nS = int(data[0][0])
    del data[0]

    for S_indx in range(nS):
            
        #new version
        if (len(data[0]) == 1 and len(data[1]) == 1):
            EF_list.append(float(data[0][0]))
            dim = int(data[1][0])
            dim_list.append(dim)
            del data[0:2]
        #old version
        else:
            dim = int(data[0][0])
            dim_list.append(dim)
            del data[0]

        Si = np.zeros((dim,dim), dtype='complex')
        for row_count in range(dim):
            row = data[:dim]
            re = np.take(row, [2], axis=1)
            im = np.take(row, [3], axis=1)
            Si[row_count] = re.flatten() + I * im.flatten()
            del data[:dim]

        S.append(Si)
     
    if (old_ver == 0):
        return np.array(S), np.array(dim_list), np.array(EF_list)
    elif (old_ver != 0):
        return np.array(S), np.array(dim_list)

    close(S_str)


def write_states(states,
                 n_states,
                 nin_Max,
                 energ,
                 filen,
                 state_sort=''):
   '''
   write states to file in proper format for Flo's code

   parameters:
   states: ndarray, float
      matrix containing row-wise states to write
   n_states: int
      number of states to be propagated by Flo's code
   nin_Max: int
      actual number of open modes
   energ: ndarray, float
      scattering energy (k**2/2)
   filen: str
      namestring of calculation
   state_sort: str (optional)
      string containing information which kind of
      states are calculated and to be propagated
   '''

   if (state_sort == ''):
      f = open('coeff.'+filen+'.dat', 'w')
   else:
      f = open('coeff.'+filen+'.'+state_sort+'.dat', 'w')

   states = scat.inflate_small_mat(states, nin_Max)
   # append zeros to make number and size of vectors match
   # nin_Max = actual number of open modes in system

   f.write('1 -1\n')
   f.write('1.0 0.5\n')
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

