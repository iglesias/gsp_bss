function wrapper_bss_nuclear

verbose_bss_nuclear = true;
verbose_self = true;
do_plot = true;

[truth, model, y] = bss_gen_problem;
[Z1_hat, Z2_hat] = bss_nuclear(y, model.A, model.G.V, verbose_bss_nuclear);

if verbose_self
  bss_print_summary(Z1_hat, Z2_hat, truth, model, y);
end
  
if do_plot
  bss_plot_results(truth.Z1, truth.Z2, Z1_hat, Z2_hat);
end
