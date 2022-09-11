function analyse_brain_pairs(predictor, probs, do_normalized_predictor)

if nargin < 3
  do_normalized_predictor = false;
end

figure
hold on

all_y = [];
all_x = [];

for i = 1:5
  x = predictor(i, (i+1):end);
  if do_normalized_predictor
    x = x / max(x);
    x = 1-x;
  end

  y = probs(i, (i+1):end);

  scatter(x, y, 'LineWidth', 2);
  
  all_x = [all_x; x'];
  all_y = [all_y; y'];
end

if do_normalized_predictor
  xlabel_str = 'Normalized predictor';
else
  xlabel_str = 'Predictor';
end

filtered_outliers_x = all_x([1:5 8:9 11:end]);
filtered_outliers_y = all_y([1:5 8:9 11:end]);

coeffs = regress(filtered_outliers_y, [ones(length(filtered_outliers_x), 1) filtered_outliers_x]);
x_axis = linspace(min(filtered_outliers_x), max(filtered_outliers_x));
plot(x_axis, coeffs(1) + coeffs(2)*x_axis, 'k--', 'LineWidth', 2)

hold off
box on
grid on
xlabel(xlabel_str, 'FontSize', 14)
ylabel('Rate of successful recovery', 'FontSize', 14)
set(gca, 'FontSize', 14)

end
