function wrapper_multigraph_bss_nuclear_direct

verbose_multigraph_bss_nuclear_direct = false;
verbose_self = true;
do_plot = true;

params.L = 3;
params.N = 50;
params.S = 1;
params.numGraphs = 2;

[truth, model, y] = multigraph_bss_gen_problem(params);
Z_hat = multigraph_bss_nuclear_direct(y, model.A, model.V, ...
                                      verbose_multigraph_bss_nuclear_direct);
if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
end

if do_plot
  plot_Zs(truth.Z, Z_hat)
end

end
