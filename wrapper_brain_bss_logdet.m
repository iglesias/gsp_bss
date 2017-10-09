function wrapper_bss_logdet_brain

verbose_multigraph_bss_logdet = true;
verbose_self = true;

[truth, model, y] = brain_bss_gen_problem([1 6]);
[Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, verbose_multigraph_bss_logdet);

if verbose_self
  multigraph_bss_print_summary(Z_hat, truth, model, y);
end

end
