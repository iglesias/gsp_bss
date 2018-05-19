fprintf('\n\n')
fprintf('Equality constraint test: %d\n', norm(y - V*A*(Z1_hat(:) + Z2_hat(:))))
fprintf('norm(Z1)=%d\n', norm(Z1))
fprintf('norm(Z2)=%d\n', norm(Z2))
fprintf('norm(Z1_hat)=%d\n', norm(Z1_hat))
fprintf('norm(Z2_hat)=%d\n', norm(Z2_hat))
fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - Z1))
fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - Z2))
fprintf('norm(Z1_hat-Z2)=%d\n', norm(Z1_hat - Z2))
fprintf('norm(Z2_hat-Z1)=%d\n', norm(Z2_hat - Z1))
fprintf('norm([Z1+Z2]-[Z1_hat+Z2_hat])=%d\n', norm([Z1+Z2] - [Z1_hat+Z2_hat], 'fro'))
fprintf('norm(Z1-Z2)=%d\n', norm(Z1 - Z2))
fprintf('norm(Z1_hat-Z2_hat)=%d\n', norm(Z1_hat - Z2_hat))
fprintf('norm_nuc(Z1)+norm_nuc(Z2)=%d\n', norm_nuc(Z1)+norm_nuc(Z2))
fprintf('norm_nuc(Z1_hat)=%d\n', norm_nuc(Z1_hat))
fprintf('norm_nuc(Z2_hat)=%d\n', norm_nuc(Z2_hat))
fprintf('norm_nuc(Z1_hat)+norm_nuc(Z2_hat)=%d\n', norm_nuc(Z1_hat)+norm_nuc(Z2_hat))
fprintf('\n\n')
fprintf('%d\n', norm(Ur(:,1)*Sr(1,1)*Vr(:,1)' - Z1))
fprintf('%d\n', norm(Ur(:,2)*Sr(2,2)*Vr(:,2)' - Z2))
fprintf('%d\n', norm(Ur(:,1)*Sr(1,1)*Vr(:,1)' - Z2))
fprintf('%d\n', norm(Ur(:,2)*Sr(2,2)*Vr(:,2)' - Z1))
fprintf('\n\n')
