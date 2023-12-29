# Analyzing the bound
#      ||conjugate(u_i,Omega) transpose(u_i,Omega) - I_N,Omega|| <= rho_U(S) (34)
# [https://arxiv.org/pdf/1604.07234.pdf]
#
# in comparison with ||a_l,Omega a_l,Omega||^2 <= max(1, ||a_l||^2_2)
# [Ling & Strohmer self-calibration]
#
# Started from github.com/iglesias/gsp_bss/blob/master/etc/blindem_linear_operators.ipynb

import random

import numpy
import numpy.matlib

N = 10
S = 3

# Sample S unique elements from 1,...,N.
Omega = random.sample(range(N), S)
s = numpy.matlib.randn(S,1)

print(f'N={N}, S={S}, Omega={Omega}')

#C_Omegat = numpy.zeros([N,S])
#C_Omegat[Omega,:] = numpy.eye(S)

#x = C_Omegat@s


import numpy as np

# Directed cycle.
GSO = np.diag(np.ones(N-1), 1)
GSO[-1,0] = 1

eigvals, V = np.linalg.eig(GSO)
V = np.matrix(V)
U = np.linalg.inv(V)
U = np.sqrt(N)*U

# project
Omegaperp = np.setdiff1d(range(N), Omega, assume_unique=False)

import matplotlib.pyplot

U_Omega = U
U_Omega[:,Omegaperp] = 0 # ojo!
matplotlib.pyplot.spy(U_Omega);

def norm2(u):
    assert(N == u.shape[0] and 1 == u.shape[1])
    # The norm in the N-dimensional complex space is
    # ||u|| := [ |u_1|^2 + |u_2|^2 + ... + |u_N|^2 ]^(1/2)
    #          square root of the sum of the squares of
    #             the moduli of the complex numbers
    #
    # Norm squared!
    return np.sum(np.array(np.abs(u))**2)

assert(N == U.shape[0])
for i in range(N):
    print('='*100)
    u_i = U[i,:].T
    print('||u_i||^2 = ', norm2(u_i), np.linalg.norm(u_i)**2)

    print('*'*100)
    A = np.conj(u_i) @ u_i.T
    print('conjugate(u_i) transpose(u_i) = ')
    matrix = A
    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix]))

print('#'*100)
b = np.linalg.norm(np.conj(U).T @ U - N*np.eye(N))
print('||U^H U - N*I_N|| =',b)

