function analysis_bss

%fname_prefix = 'bss_random_filters_*';
%fname_prefix = 'bss_hLP01_hHP-10_filters_*';
%fname_prefix = 'bss_hLP10_hHP01_filters_*';
fname_prefix = 'bss_nuclear_hLP10_hHP01_filters_*';

files = dir(fname_prefix);
eq_constraint_test = zeros(length(files), 1);
Zi_err = zeros(length(files), 2);
Z_err = zeros(length(files), 1);
Zi_hat_rank = zeros(length(files), 2);
Zi_diff = zeros(length(files), 2);
for i = 1:length(files)
    load(files(i).name)
    eq_constraint_test(i) = norm(y - V*A*(Z1_hat(:) + Z2_hat(:)));
    Zi_err(i, :) = [norm(Z1-Z1_hat, 'fro') norm(Z2-Z2_hat, 'fro')];
    Z_err(i) = norm((Z1+Z2) - (Z1_hat+Z2_hat), 'fro');
    Zi_hat_rank(i, :) = [rank(Z1_hat) rank(Z2_hat)];
    Zi_diff(i, :) = [norm(Z1-Z2, 'fro') norm(Z1_hat-Z2_hat, 'fro')];
end

figure(1), clf
stem(eq_constraint_test)
minmax_values = minmax(eq_constraint_test);
title(sprintf('Constraint test \nmin=%d max=%d', minmax_values(1), ...
              minmax_values(2)))
xlabel('Simulation index')

figure(2), clf
hold on
plot(Zi_err(:, 1), 'bo-', 'LineWidth', 2)
plot(Zi_err(:, 2), 'mx-', 'LineWidth', 2)
plot(Z_err, 'ks-.')
hold off
title('Z estimation errors')
legend('norm(Z1-Z1\_hat)', 'norm(Z2-Z2\_hat)', 'norm(Z-Z\_hat)')
xlabel('Simulation index')

if all(Zi_hat_rank(:) == 2)
    fprintf('All Z_hat have rank equal to 2.\n')
else
    warning('Rank of Z_hat different to 2!\n')
end

figure(3), clf
hold on
plot(Zi_diff(:, 1), 'bo-.')
plot(Zi_diff(:, 2), 'mx-.')
hold off
title('Difference between Z1 and Z2 and Z1\_hat and Z2\_hat')
legend('norm(Z1-Z2)', 'norm(Z1\_hat-Z2\_hat)')
xlabel('Simulation index')

figure(4), clf
plot(abs(Zi_err(:, 1)-Zi_err(:,2)), 'bo-.')
title('Difference between Zi estimation errors')
legend('abs(norm(Z1-Z1\_hat) - norm(Z2-Z2\_hat))')
xlabel('Simulation index')

end
