function wrapper_multigraph_bss_logdet

verbose_multigraph_bss_logdet = false;
verbose_self = true;
do_plot = false;

params.L = 5;
params.N = 34;
params.S = 1;
params.numGraphs = 2;
fprintf('L=%d N=%d S=%d\n', params.L, params.N, params.S);

[truth, model, y] = multigraph_bss_gen_problem(params);
Z_hat = multigraph_bss_logdet(y, model.A, model.V, ...
                              verbose_multigraph_bss_logdet);
if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
end

if do_plot
  plot_Zs(truth.Z, Z_hat)
end

end
