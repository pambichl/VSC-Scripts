pht = np.array([[
            phasetimes.calc_phase_times(S[i,j], dims, dk)
            for j in range(nconf)] for i in range(nr)])
pht_means = np.array([[
            phasetimes.calc_phase_means(pht[i,j], dims)
            for j in range(nconf)] for i in range(nr)])
means_phase = np.mean(pht_means, axis=1)
            

#print pht[0,0]
print means_phase / means



meansIm = np.array([[np.mean(q_mean[i,:,j].imag) for j in range(nfreq)] for i in range(nr)])
stddevsIm = np.array([[np.std(q_mean[i,:,j].imag)  for j in range(nfreq)] for i in range(nr)])
print "number of non-unitary Q matrices: %i"%herm_counter



        #
        #print ut.sort_eig(scat.calc_Q(s,dk)[0])[0]
        #
        #q_mean = np.trace(scat.calc_Q(s,dk)[0] +
        #                  1J/2. * (np.dot(s[1].T.conj(), k_mat[w]) - np.dot(k_mat[w], s[1]))) / (2*dims[w])



# for w in range(3*nfreq):
#     energ = energs[w]
#     energs_y = 1./dy**2 * (1.-np.cos(np.arange(1,dims[int(w/3)]+1,1)*np.pi/W*dy))
#     energs_x = energ - energs_y
#     k_x = np.arccos(1.-energs_x*dx**2)/dx
#     k_x_list.append(np.array(k_x))
#     #k_x = np.sqrt(2*energ-(np.arange(1,dims[w]+1,1)*np.pi/W*dy)**2)
#     k_block = np.diag(1./k_x**2)
#     zero = np.zeros((dims[int(w/3)],dims[int(w/3)]), dtype="complex")
#     k_mat.append(np.bmat([[k_block, zero], [zero, k_block]]))

#print "kx:", k_x_list
#print k_mat[0]
#print np.diag(k_mat[0])
#print np.sqrt(2*energs[0]-(np.arange(1,dims[w]+1,1)*np.pi/W*dy)**2)



#    for prob in problems:
#        subprocess.call(["rm", data_direc+"Smat."+prob+".dat"])



W = 10.0
dy= W / 150.0
dx = dy
k_x_list = []
k_mat = []
