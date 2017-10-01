function wrapper_bss_logdet_freq_svd

verbose_bss_logdet_freq = true;
verbose_self = true;
do_plot = false;

[truth, model, y_tilde] = bss_gen_freq_problem;
[Z_hat, iter] = bss_logdet_freq(y_tilde, model.A, verbose_bss_logdet_freq);

[UZPsi, SZPsi, VZPsi] = svd(Z_hat);
ZPsi1_hat = SZPsi(1,1) * UZPsi(:,1) * VZPsi(:,1)';
ZPsi2_hat = SZPsi(2,2) * UZPsi(:,2) * VZPsi(:,2)';

if verbose_self
  bss_print_freq_summary(ZPsi1_hat, ZPsi2_hat, truth, model);
end

if do_plot
  bss_plot_results(truth.ZPsi1, truth.ZPsi2, ZPsi1_hat, ZPsi2_hat);
end

end
