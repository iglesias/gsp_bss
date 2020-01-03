%SS = [1 2 3];
%NUM_FILTERS = [2 3];
%NOISE = [1e-6 1e-5 1e-4 1e-3 1e-2];

SS = [1];
NUM_FILTERS = 2;
NOISE = [1e-6 1e-4 1e-2];

num_plots = length(SS)*length(NUM_FILTERS);
plot_mean_data = zeros(length(NOISE), num_plots);
plot_median_data = zeros(length(NOISE), num_plots);
fname_pattern = 'singlegraph_karate_S%d_numFilters%dnoise%d.mat';

for n = 1:length(NOISE)
  plot_idx = 1;
  for S = SS
    for numFilters = NUM_FILTERS
      fname = sprintf(fname_pattern, S, numFilters, NOISE(n));
      fprintf('Reading %s\n', fname);
      load(fname)
      SS = [1];
      NUM_FILTERS = 2;
      NOISE = [1e-6 1e-4 1e-2];

      plot_mean_data(n, plot_idx) = mean(recovery_performance);
      plot_median_data(n, plot_idx) = median(recovery_performance);

      plot_idx = plot_idx + 1;
    end
  end
end

figure
semilogy(NOISE, plot_median_data, 's-')
%legend('F = 2, S = 3', 'F = 3, S = 3', 'F = 2, S = 6', 'F = 3, S = 6', 'Location', 'NorthEast')
xlabel('\sigma_n')
ylabel('Recovery performance (median)')
grid on
title('Karate network')

figure
semilogy(NOISE, plot_mean_data, 's-')
%legend('F = 2, S = 3', 'F = 3, S = 3', 'F = 2, S = 6', 'F = 3, S = 6', 'Location', 'NorthEast')
xlabel('\sigma_n')
ylabel('Recovery performance (mean)')
grid on
title('Karate network')
