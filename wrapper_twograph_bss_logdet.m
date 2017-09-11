clearvars, twograph_bss_gen_problem
verbose_twograph_bss_logdet = true;
verbose_self = true;
[Z1_hat, Z2_hat] = twograph_bss_logdet(y, A1, G(1).V, A2, G(2).V, verbose_twograph_bss_logdet);
if verbose_self, twograph_bss_print_summary, end
