from math import sqrt
from numpy.linalg import norm

def recovery_assessment(Z, Z_hat):
  assert(len(Z) == Z_hat.shape[2])

  num = 0
  den = 0

  for p in range(len(Z)):
    num += norm(Z_hat[:,:,p] - Z[p])**2
    den += norm(Z[p])**2

  return sqrt(num) / sqrt(den)

if __name__ == '__main__':
  from loadmat import loadmat
  mat = loadmat('20180124/singlegraph_bss_logdet_knownx_02bJ6ZgjQVAdrHUsO6xJ.mat')
  print(recovery_assessment(mat['truth']['Z'], mat['Z_hat']))
