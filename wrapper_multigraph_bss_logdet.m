function wrapper_multigraph_bss_logdet

verbose_multigraph_bss_logdet = false;
verbose_self = false;

params.numGraphs = 5;
params.L = 3;
params.N = 100;
%params.S = 5;

[truth, model, y] = multigraph_bss_gen_problem(params);
Z_hat = multigraph_bss_logdet(y, model.A, model.V, ...
                              verbose_multigraph_bss_logdet);
if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
end

end
