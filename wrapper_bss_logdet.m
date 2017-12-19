function wrapper_bss_logdet

verbose_bss_logdet = false;
verbose_self = true;
do_plot = false;

[truth, model, y] = singlegraph_svd_bss_gen_problem;
[Z1_hat, Z2_hat] = bss_logdet(y, model.A, model.G.V, verbose_bss_logdet);

if verbose_self
  bss_print_summary(Z1_hat, Z2_hat, truth, model, y);
end

if do_plot
  bss_plot_results(truth.Z{1}, truth.Z{2}, Z1_hat, Z2_hat);
end

end
