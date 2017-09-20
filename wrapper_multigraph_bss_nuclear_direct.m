function wrapper_multigraph_bss_nuclear_direct

verbose_multigraph_bss_nuclear_direct = true;
verbose_self = false;

[truth, model, y] = multigraph_bss_gen_problem;
Z_hat = multigraph_bss_nuclear_direct(y, model.A, model.V, ...
                                      verbose_multigraph_bss_nuclear_direct);
if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
end

end
