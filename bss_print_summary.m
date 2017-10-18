function bss_print_summary(Z1_hat, Z2_hat, truth, model, y)
%BSS_PRINT_SUMMARY Suitable only for the single-graph, two-filter problem case.

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', ...
        norm(y - model.G.V*model.A*(Z1_hat(:) + Z2_hat(:))))
fprintf('norm(Z1)=%d\n', norm(truth.Z{1}))
fprintf('norm(Z2)=%d\n', norm(truth.Z{2}))
fprintf('norm(Z1_hat)=%d\n', norm(Z1_hat))
fprintf('norm(Z2_hat)=%d\n', norm(Z2_hat))
fprintf('Recovery assessment: %d\n', ...
        recovery_assessment({truth.Z{1}, truth.Z{2}}, cat(3, Z1_hat, Z2_hat)))
fprintf('Recovery assessment: %d\n', ...
        recovery_assessment({truth.Z{1}, truth.Z{2}}, cat(3, Z2_hat, Z1_hat)))
fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - truth.Z{1}))
fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - truth.Z{2}))
fprintf('norm(Z1_hat-Z2)=%d\n', norm(Z1_hat - truth.Z{2}))
fprintf('norm(Z2_hat-Z1)=%d\n', norm(Z2_hat - truth.Z{1}))
fprintf('norm([Z1+Z2]-[Z1_hat+Z2_hat])=%d\n', ...
        norm([truth.Z{1}+truth.Z{2}] - [Z1_hat+Z2_hat], 'fro'))
fprintf('norm(Z1-Z2)=%d\n', norm(truth.Z{1} - truth.Z{2}))
fprintf('norm(Z1_hat-Z2_hat)=%d\n', norm(Z1_hat - Z2_hat))
fprintf('norm_nuc(Z1)+norm_nuc(Z2)=%d\n', norm_nuc(truth.Z{1})+norm_nuc(truth.Z{2}))
fprintf('norm_nuc(Z1_hat)=%d\n', norm_nuc(Z1_hat))
fprintf('norm_nuc(Z2_hat)=%d\n', norm_nuc(Z2_hat))
fprintf('norm_nuc(Z1_hat)+norm_nuc(Z2_hat)=%d\n', ...
        norm_nuc(Z1_hat)+norm_nuc(Z2_hat))
fprintf('rank(Z1)=%d rank(Z2)=%d\n', rank(truth.Z{1}), rank(truth.Z{2}))
fprintf('\n\n')

end
