function analysis_alpha

alphas = 0:0.1:1.0;

Zi_hat_norms = zeros(length(alphas), 2);
Zi_hat_diff = zeros(length(alphas), 1);

for i = 1:length(alphas)
  data = analysis_helper('alpha', alphas(i));
  Zi_hat_norms(i, :) = mean(data.Zi_hat_norms, 1);
  Zi_hat_diff(i) = mean(data.Zi_hat_diff);
end

plot(alphas, Zi_hat_norms, 'o-', alphas, Zi_hat_diff, 'x--')
xlabel('\alpha')
ylabel('||A||_F')
legend('||Z1\_hat||_F', '||Z2\_hat||_F', '||Z1\_hat-Z2\_hat||_F', 'Location', 'best')
title(sprintf('Recovered matrices with objective containing \n \\alpha ||Z1||_* + (1-\\alpha) ||Z2||_*'))
grid on

end
