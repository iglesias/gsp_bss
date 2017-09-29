function bss_print_summary(Z1_hat, Z2_hat, truth, model, y)

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', ...
        norm(y - model.G.V*model.A*(Z1_hat(:) + Z2_hat(:))))
fprintf('norm(Z1)=%d\n', norm(truth.Z1))
fprintf('norm(Z2)=%d\n', norm(truth.Z2))
fprintf('norm(Z1_hat)=%d\n', norm(Z1_hat))
fprintf('norm(Z2_hat)=%d\n', norm(Z2_hat))
fprintf('Recovery assessment: %d\n', ...
        recovery_assessment({truth.Z1, truth.Z2}, cat(3, Z1_hat, Z2_hat)))
fprintf('Recovery assessment: %d\n', ...
        recovery_assessment({truth.Z1, truth.Z2}, cat(3, Z2_hat, Z1_hat)))
fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - truth.Z1))
fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - truth.Z2))
fprintf('norm(Z1_hat-Z2)=%d\n', norm(Z1_hat - truth.Z2))
fprintf('norm(Z2_hat-Z1)=%d\n', norm(Z2_hat - truth.Z1))
fprintf('norm([Z1+Z2]-[Z1_hat+Z2_hat])=%d\n', ...
        norm([truth.Z1+truth.Z2] - [Z1_hat+Z2_hat], 'fro'))
fprintf('norm(Z1-Z2)=%d\n', norm(truth.Z1 - truth.Z2))
fprintf('norm(Z1_hat-Z2_hat)=%d\n', norm(Z1_hat - Z2_hat))
fprintf('norm_nuc(Z1)+norm_nuc(Z2)=%d\n', norm_nuc(truth.Z1)+norm_nuc(truth.Z2))
fprintf('norm_nuc(Z1_hat)=%d\n', norm_nuc(Z1_hat))
fprintf('norm_nuc(Z2_hat)=%d\n', norm_nuc(Z2_hat))
fprintf('norm_nuc(Z1_hat)+norm_nuc(Z2_hat)=%d\n', ...
        norm_nuc(Z1_hat)+norm_nuc(Z2_hat))
fprintf('rank(Z1)=%d rank(Z2)=%d\n', rank(truth.Z1), rank(truth.Z2))
fprintf('svd(Z1_hat)=(%d, %d)\n', svd(Z1_hat))
fprintf('svd(Z2_hat)=(%d, %d)\n', svd(Z2_hat))
fprintf('\n\n')

end
