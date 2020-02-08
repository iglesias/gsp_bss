clearvars

SS = flip([1 2 3]);
NUM_FILTERS = flip([2 3]);
NOISE = [1e-6 1e-5 1e-4 1e-3 1e-2];

num_plots = length(SS)*length(NUM_FILTERS);
plot_mean_data = zeros(length(NOISE), num_plots);
plot_median_data = zeros(length(NOISE), num_plots);
fname_pattern = 'singlegraph_karate_S%d_numFilters%dnoise%d';

for n = 1:length(NOISE)
  plot_idx = 1;
  for numFilters = NUM_FILTERS
    for S = SS
      fname = sprintf(fname_pattern, S, numFilters, NOISE(n));
      fprintf('Reading %s\n', fname);
      load(fname, '-mat', 'recovery_performance')

%       plot_mean_data(n, plot_idx) = mean(recovery_performance);
      plot_median_data(n, plot_idx) = median(recovery_performance);

      plot_idx = plot_idx + 1;
      if (n == 1)
        if ~exist('legend_str', 'var')
          legend_str = {sprintf('F = %d, S = %d', numFilters, S)};
        else
          legend_str = [legend_str, {sprintf('F = %d, S = %d', numFilters, S)}];
        end
      end
    end
  end
end

figure
loglog(NOISE, plot_median_data, 's-', 'LineWidth', 1.5)
% legend('F = 2, S = 1', 'F = 3, S = 1', 'F = 2, S = 2', 'F = 3, S = 2', 'F = 2, S = 3', 'F = 3, S = 3', 'Location', 'NorthEastOutside')
legend(legend_str, 'Location', 'NorthEastOutside')
xlabel('\sigma_n')
ylabel('Recovery performance (median)')
grid on
title('Karate network')

% figure
% loglog(NOISE, plot_mean_data, 's-', 'LineWidth', 1.5)
% legend('F = 2, S = 3', 'F = 3, S = 3', 'F = 2, S = 6', 'F = 3, S = 6', 'F = 2, S = 3', 'F = 3, S = 3', 'Location', 'NorthEastOutside')
% xlabel('\sigma_n')
% ylabel('Recovery performance (mean)')
% grid on
% title('Karate network')
