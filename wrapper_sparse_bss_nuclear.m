verbose_sparse_bss_nuclear = false;
verbose_self = false;
taux = 0.5;
[Z1_hat, Z2_hat] = sparse_bss_nuclear(y, A, V, [], [], verbose_sparse_bss_nuclear);
if verbose_self, bss_print_summary, end
