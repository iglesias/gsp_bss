function wrapper_brain_bss_logdet

verbose_multigraph_bss_logdet = true;
verbose_self = true;

params.numGraphs = 2;
params.L = 2;
params.S = 1;
params.brain_idxs = [2 3];

[truth, model, y] = brain_bss_gen_problem(params);
Z_hat = multigraph_bss_logdet(y, model.A, model.V, ...
                              verbose_multigraph_bss_logdet);

if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
  plot_Zs(truth.Z, Z_hat);
end

end
