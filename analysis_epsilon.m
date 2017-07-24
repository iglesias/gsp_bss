function analysis_epsilon

%epsilons = 0:0.1:1.0;
epsilons = [1e-10 1e-9 1e-8 1e-7 1e-4];
Zi_hat_norms = zeros(length(epsilons), 2);
Zi_hat_diff = zeros(length(epsilons), 1);

for i = 1:length(epsilons)
  data = analysis_helper('epsilon', epsilons(i));
  eq_constraint_test(i, :) = data.eq_constraint_test;
  Zsum_test(i, :) = data.Zsum_test;
  Zi_hat_norms(i, :) = mean(data.Zi_hat_norms, 1);
  Zi_hat_diff(i) = mean(data.Zi_hat_diff);
  Z1_col_norms(i, :) = mean(data.Z1_col_norms, 1);
  Z1_hat_col_norms(i, :) = mean(data.Z1_hat_col_norms, 1);
  Z1_hat_col_norms_std(i, :) = std(data.Z1_hat_col_norms);
  Z2_col_norms(i, :) = mean(data.Z2_col_norms, 1);
  Z2_hat_col_norms(i, :) = mean(data.Z2_hat_col_norms, 1);
  Z2_hat_col_norms_std(i, :) = std(data.Z2_hat_col_norms);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

%plot(epsilons, Zi_hat_norms, 'o-')
%legend('||Z1\_hat||_F', '||Z2\_hat||_F', 'Location', 'best')

semilogx(epsilons, [Zi_hat_norms Zi_hat_diff], 'o-')
legend('||Z1\_hat||_F', '||Z2\_hat||_F', '||Z1\_hat-Z2\_hat||_F', 'Location', 'best')

xlabel('\epsilon')
ylabel('||A||_F')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

subplot(221)
semilogx(epsilons, Z1_col_norms, 'o-')
ylim([0 2])
title('Z1')
grid on

subplot(223)
semilogx(epsilons, Z1_hat_col_norms, 'o-', ...
         epsilons, [Z1_hat_col_norms+Z1_hat_col_norms_std Z1_hat_col_norms-Z1_hat_col_norms_std], 'k-.')
ylim([0 2])
title('Z1\_hat')
grid on

subplot(222)
semilogx(epsilons, Z2_col_norms, 'o-')
ylim([0 2])
title('Z2')
grid on

subplot(224)
semilogx(epsilons, Z2_hat_col_norms, 'o-', ...
         epsilons, [Z2_hat_col_norms+Z2_hat_col_norms_std Z2_hat_col_norms-Z2_hat_col_norms_std], 'k-.')
ylim([0 2])
title('Z2\_hat')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
boxplot(Zsum_test', 'Labels', epsilons)
xlabel('\epsilon')
ylabel('||(Z1+Z2)-(Z1\_hat+Z2\_hat)||_F')

end
