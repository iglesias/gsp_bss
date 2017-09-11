function twograph_bss_print_summary(Z1_hat, Z2_hat, truth, model, y)

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', norm(y - (model.G(1).V*model.A1*Z1_hat(:) + model.G(2).V*model.A2*Z2_hat(:))))
fprintf('Recovery assessment: %d\n', recovery_assessment(truth.Z1, truth.Z2, Z1_hat, Z2_hat))
fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - truth.Z1))
fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - truth.Z2))
fprintf('rank(Z1)=%d rank(Z2)=%d\n', rank(truth.Z1), rank(truth.Z2))
fprintf('svd(Z1_hat)=(%d, %d)\n', svd(Z1_hat))
fprintf('svd(Z2_hat)=(%d, %d)\n', svd(Z2_hat))
fprintf('\n\n')

end
