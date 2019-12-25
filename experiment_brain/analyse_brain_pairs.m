function analyse_brain_pairs(predictor, probs)

for i = 1:1

x = predictor(i, (i+1):end);
x = x / max(x);
x = 1-x;

figure
[probs_i, idxs] = sort(probs(i, (i+1):end));
plot(x(idxs), probs_i, 'o--', 'LineWidth', 2);
box on
grid on
xlabel('Normalized predictor')
ylabel('Success probability')

end

end
