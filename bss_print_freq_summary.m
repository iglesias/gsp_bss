function bss_print_freq_summary(ZPsi1_hat, ZPsi2_hat, truth, model)

assert(all(size(ZPsi1_hat) == size(ZPsi2_hat)))
svd_fmt_str = prepare_svd_fmt_str(size(ZPsi1_hat, 2));

fprintf('\n\n')
fprintf('Recovery assessment: %d\n', ...
        recovery_assessment({truth.ZPsi1, truth.ZPsi2}, cat(3, ZPsi1_hat, ZPsi2_hat)))
fprintf('Recovery assessment: %d\n', ...
        recovery_assessment({truth.ZPsi1, truth.ZPsi2}, cat(3, ZPsi2_hat, ZPsi1_hat)))
fprintf('norm([ZPsi1+ZPsi2]-[ZPsi1_hat+ZPsi2_hat])=%d\n', ...
        norm([truth.ZPsi1+truth.ZPsi2] - [ZPsi1_hat+ZPsi2_hat], 'fro'))
fprintf('rank(ZPsi1)=%d rank(ZPsi2)=%d\n', rank(truth.ZPsi1), rank(truth.ZPsi2))
fprintf(sprintf('svd(ZPsi1_hat)=(%s)\n', svd_fmt_str), svd(ZPsi1_hat))
fprintf(sprintf('svd(ZPsi2_hat)=(%s)\n', svd_fmt_str), svd(ZPsi2_hat))
fprintf('\n\n')

end
