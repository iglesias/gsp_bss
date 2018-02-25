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

num_success = 0

mat_dir = '20180125'
mat_files = [f for f in listdir(mat_dir) if isfile(join(mat_dir, f))]

for mat_file in mat_files:
  #print(mat_file)
  mat = loadmat(join(mat_dir, mat_file))

  try:
    N, L, P = mat['Z_hat'].shape
  except ValueError:
    continue

  S = len(flatnonzero(mat['truth']['x'][:,0]))

  if P != 2: continue

  for p in range(P):
    if isinstance(mat['nonzero_idxs'][p], int):
      mat['nonzero_idxs'][p] = array([mat['nonzero_idxs'][p]])

  zero_known_len = len(mat['zero_idxs'][0])
  nonzero_known_len = len(mat['nonzero_idxs'][0])
  known_ratio = (zero_known_len + nonzero_known_len) / N

  perf = recovery_assessment(mat['truth']['Z'], mat['Z_hat'])

  #TODO probably needs tweaking if P > 2.
  assert(P == 2)
  zero_intersection_len = len(set(mat['zero_idxs'][0]).intersection(mat['zero_idxs'][1]))
  nonzero_intersection_len = len(set(mat['nonzero_idxs'][0]).intersection(mat['nonzero_idxs'][1]))
  assert(P == 2)

  if perf < 1e-3:
    zero_intersection_hist.add(zero_intersection_len)
    nonzero_intersection_hist.add(nonzero_intersection_len)
    S_hist.add(S)
    num_success += 1

  print('\t', 'K=%.2f' % known_ratio, 'L=%d' % L, 'N=%d' % N, 'P=%d' % P, 'S=%d' % S, perf)

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
