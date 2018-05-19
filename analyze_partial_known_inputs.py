# System imports.
from os import listdir
from os.path import isfile, join

# Other libraries imports.
from numpy import array, flatnonzero

# Own imports.
from loadmat import loadmat
from recovery_assessment import recovery_assessment

class IntegerHistogram:
  def __init__(self):
    self._dict = {}

  def add(self, n):
    if n in self._dict.keys():
      self._dict[n] += 1
    else:
      self._dict[n] = 1

zero_intersection_hist = IntegerHistogram()
nonzero_intersection_hist = IntegerHistogram()
S_hist = IntegerHistogram()

zero_intersection_success_hist = IntegerHistogram()
nonzero_intersection_success_hist = IntegerHistogram()
S_success_hist = IntegerHistogram()

num_success = 0

mat_dir = '20180124'
mat_files = [f for f in listdir(mat_dir) if isfile(join(mat_dir, f))]

nonzero_known_success = 0

for mat_file in mat_files:
  #print(mat_file)
  mat = loadmat(join(mat_dir, mat_file))

  # Obtain #nodes (N), #filter coefficients (L), and #filters (P), handling the P=1 case.
  N, L = mat['Z_hat'].shape[0:2]
  P = mat['Z_hat'].shape[2] if len(mat['Z_hat'].shape) > 2 else 1

  if P != 2: continue

  S = len(flatnonzero(mat['truth']['x'][:,0]))

  for p in range(P):
    if isinstance(mat['nonzero_idxs'][p], int):
      mat['nonzero_idxs'][p] = array([mat['nonzero_idxs'][p]])

  zero_known_len = len(mat['zero_idxs'][0])
  nonzero_known_len = len(mat['nonzero_idxs'][0])
  known_ratio = 1.0*(zero_known_len + nonzero_known_len) / N

  perf = recovery_assessment(mat['truth']['Z'], mat['Z_hat'])

  #TODO probably needs tweaking if P > 2.
  # Think about the meaning of intersection with P>2:
  # - intersection of all sources?
  # - pair-wise intersections?
  assert(P == 2)
  zero_intersection_len = len(set(mat['zero_idxs'][0]).intersection(mat['zero_idxs'][1]))
  nonzero_intersection_len = len(set(mat['nonzero_idxs'][0]).intersection(mat['nonzero_idxs'][1]))
  assert(P == 2)

  S_hist.add(S)
  zero_intersection_hist.add(zero_intersection_len)
  nonzero_intersection_hist.add(nonzero_intersection_len)

  if perf < 1e-3:
    if len(mat['nonzero_idxs'][0]) + len(mat['nonzero_idxs'][1]) >= 1:
      nonzero_known_success += 1

    zero_intersection_success_hist.add(zero_intersection_len)
    nonzero_intersection_success_hist.add(nonzero_intersection_len)
    S_success_hist.add(S)
    num_success += 1

#  print('\t K=%.2f L=%d N=%d P=%d S=%d: %.0e' % (known_ratio, L, N, P, S, perf))

#print(mat.keys())
#
#print(mat['nonzero_idxs'])
#print(mat['zero_idxs'])
#
##mat['model']
#print(mat['model'].keys())
#
#print(mat['iter'])
#
#print(mat['truth']['h'])
#print(mat['truth']['x'])
