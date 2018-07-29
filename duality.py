import numpy as np

def make_connected_erdos_rendi(N=10, p=0.1):
  while True:
    S = np.array(np.random.rand(N, N) < p, dtype=float)
    S = np.triu(S, k=1)
    S = S + S.T
    # Make sure it is connected.
    L = np.diag(np.sum(S, 0)) - S
    eigvals = np.linalg.eig(L)[0]
    if sum(abs(eigvals) < 1e-6) == 1:
      return S

def make_filter_matrix(S, h):
  H = h[0]*np.eye(*S.shape)
  for l in range(1, len(h)):
    H += h[l]*np.linalg.matrix_power(S, l)
  return H

N = 10
L = 4

h = np.random.randn(L, 1)
x = np.random.randn(N, 1)

S = make_connected_erdos_rendi()

eigvals, V = np.linalg.eig(S)
V = np.matrix(V)
U = np.linalg.inv(V)
Psi = np.vander(eigvals, len(h), increasing=True)
Psi = np.matrix(Psi)

H = make_filter_matrix(S, h)

# Alternative definition of H.
Htenate = V*np.matrix(np.diag(np.squeeze(np.array(Psi*h))))*U
d1 = np.linalg.norm(H - Htenate)
if d1 > 1e-10:
  print(d1)

a = V*np.matrix(np.diag(np.squeeze(np.array(U*x))))*Psi
b = Htenate
d2 = np.linalg.norm(a*h - b*x)
if d2 > 1e-10:
  print(d2)
