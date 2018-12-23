function wrapper_singlegraph_bss_logdet_knownx

verbose_optimizer = false;
verbose_self = false;

params.verbose = false;
params.numFilters = 2;
params.L = 3;
params.N = 50;
params.S = 3;

[truth, model, y] = singlegraph_bss_gen_problem(params);

known_ratio = 1;
[Z_hat, iter] = singlegraph_bss_logdet_knownx(truth.x, y, model.A, model.G.V, ...
                                              known_ratio, verbose_optimizer);

if verbose_self
  plot_Zs(truth.Z, Z_hat)
  singlegraph_bss_print_summary(Z_hat, truth, model, y);
end

end
