function analysis_epsilon

%epsilons = 0:0.1:1.0;
epsilons = [1e-10 1e-9 1e-8 1e-7 1e-4];
Zi_hat_norms = zeros(length(epsilons), 2);
Zi_hat_diff = zeros(length(epsilons), 1);

for i = 1:length(epsilons)
  data = analysis_helper('epsilon', epsilons(i));
  Zi_hat_norms(i, :) = mean(data.Zi_hat_norms, 1);
  Zi_hat_diff(i) = mean(data.Zi_hat_diff);
end

%plot(epsilons, Zi_hat_norms, 'o-')
semilogx(epsilons, [Zi_hat_norms Zi_hat_diff], 'o-')
xlabel('\epsilon')
ylabel('||A||_F')
legend('||Z1\_hat||_F', '||Z2\_hat||_F', '||Z1\_hat-Z2\_hat||_F', 'Location', 'best')
title('TODO')
grid on

end
