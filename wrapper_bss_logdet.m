verbose_sparse_bss_logdet = true;
verbose_self = true;
[Z1_hat, Z2_hat] = sparse_bss_logdet(y, A, V, [], [], verbose_sparse_bss_logdet);
if verbose_self, bss_print_summary, bss_plot_results, end
