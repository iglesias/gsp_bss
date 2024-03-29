function wrapper_bss_nuclear

verbose_bss_nuclear = false;
verbose_self = false;
do_plot = false;

params.numFilters = 2;

[truth, model, y] = singlegraph_svd_bss_gen_problem(params);
[Z1_hat, Z2_hat] = bss_nuclear(y, model.A, model.G.V, verbose_bss_nuclear);
Z_hat = cat(3, Z1_hat, Z2_hat);

if verbose_self
  do_perms = true;
  singlegraph_bss_print_summary(Z_hat, truth, model, y, do_perms);
end
  
if do_plot
  plot_Zs(truth.Z, Z_hat)
end
