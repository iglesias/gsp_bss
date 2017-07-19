function analysis_epsilon

epsilons = 0:0.1:1.0;
Zi_hat_norms = zeros(length(epsilons), 2);

for i = 1:length(epsilons)
  data = analysis_helper('epsilon', epsilons(i));
  Zi_hat_norms(i, :) = mean(data.Zi_hat_norms, 1);
end

plot(epsilons, Zi_hat_norms, 'o-')
xlabel('\epsilon')
ylabel('||Zi\_hat||_F')
legend('||Z1\_hat||_F', '||Z2\_hat||_F', 'Location', 'best')
title('TODO')
grid on

end
