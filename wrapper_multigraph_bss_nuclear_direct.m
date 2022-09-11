function wrapper_multigraph_bss_nuclear_direct

verbose_multigraph_bss_nuclear_direct = false;
verbose_self = true;

params.L = 3;
params.N = 200;
params.S = 1;
params.numGraphs = 2;
fprintf('L=%d N=%d S=%d\n', params.L, params.N, params.S);

[truth, model, y] = multigraph_bss_gen_problem;
Z_hat = multigraph_bss_nuclear_direct(y, model.A, model.V, ...
                                      verbose_multigraph_bss_nuclear_direct);
if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
end

end
