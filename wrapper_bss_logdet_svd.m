function wrapper_bss_logdet_svd

verbose_bss_logdet = false;
verbose_self = true;
do_plot = false;

[truth, model, y] = bss_gen_problem;
[Z_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);

[UZ, SZ, VZ] = svd(Z_hat, 'econ');
Z1_hat = SZ(1,1)*UZ(:,1)*VZ(:,1)';
Z2_hat = SZ(2,2)*UZ(:,2)*VZ(:,2)';

if verbose_self
  bss_print_summary(Z1_hat, Z2_hat, truth, model, y);
end

if do_plot
  bss_plot_results(truth.Z1, truth.Z2, Z1_hat, Z2_hat);
end

end
