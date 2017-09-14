function twograph_bss_print_summary(Z1_hat, Z2_hat, truth, model, y)

assert(all(size(Z1_hat) == size(Z2_hat)))
svd_fmt_str = prepare_svd_fmt_str(size(Z1_hat, 2));
y_hat = model.G(1).V*model.A1*Z1_hat(:) + model.G(2).V*model.A2*Z2_hat(:);

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', norm(y - y_hat))
fprintf('Recovery assessment: %d\n', recovery_assessment(truth.Z1, truth.Z2, Z1_hat, Z2_hat))
fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - truth.Z1))
fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - truth.Z2))
fprintf('rank(Z1)=%d rank(Z2)=%d\n', rank(truth.Z1), rank(truth.Z2))
fprintf(sprintf('svd(Z1_hat)=(%s)\n', svd_fmt_str), svd(Z1_hat))
fprintf(sprintf('svd(Z2_hat)=(%s)\n', svd_fmt_str), svd(Z2_hat))
fprintf('\n\n')

end
