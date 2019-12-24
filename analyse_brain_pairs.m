function analyse_brain_pairs(predictor, probs)

for i = 1:4

x = predictor(i, (i+1):end);
x = x / max(x);
x = 1-x;

figure
hold on
plot(i+1:6, x, 'o--', 'LineWidth', 2)
plot(i+1:6, probs(i, (i+1):end), 'o--', 'LineWidth', 2);
hold off
box on
grid on
title(sprintf('i=%d', i))
legend('Predictor', 'Success prob.', 'Location', 'SouthEast')
xlabel('j indices')
ylabel('Probability / predictor')

end

end