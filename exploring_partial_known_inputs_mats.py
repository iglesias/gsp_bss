# coding: utf-8
loadmat
from loadmat import loadmat
loadmat
dir(loadmat)
matfile = loadmat('20180124/singlegraph_bss_logdet_knownx_02bJ6ZgjQVAdrHUsO6xJ.mat')
get_ipython().magic('who ')
get_ipython().magic('whos ')
matfile
matfile.keys()
matfile['nonzero_idxs']
matfile['zero_idxs']
matfile['model']
matfile['model'].keys()
matfile.keys()
matfile['iter']
matfile['Z_hat']
matfile['Z_hat'].shape
matfile['truth']['h']
matfile['truth']['x']
