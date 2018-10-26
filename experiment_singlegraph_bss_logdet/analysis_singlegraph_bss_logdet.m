NN = [20 40 60 80 100 120];
SS = [3 6];
NUM_FILTERS = [2 3];

num_plots = length(SS)*length(NUM_FILTERS);
plot_mean_data = zeros(length(NN), num_plots);
plot_median_data = zeros(length(NN), num_plots);
fname_pattern = 'play_singlegraph_bss_logdet_N%d_S%d_L3_numFilters%d';

for n = 1:length(NN)
  plot_idx = 1;
  for S = SS
    for numFilters = NUM_FILTERS
      fname = sprintf(fname_pattern, NN(n), S, numFilters);
      fprintf('Reading %s\n', fname);
      load(fname)

      if NN(n) < 80
        %figure
        %hist(recovery_performance)
        title(sprintf('N = %d F = %d S = %d', NN(n), numFilters, S))
        xlabel('Recovery performance')
      end

      plot_mean_data(n, plot_idx) = mean(recovery_performance);
      plot_median_data(n, plot_idx) = median(recovery_performance);

      plot_idx = plot_idx + 1;
    end
  end
end

figure
semilogy(NN, plot_median_data, 's-')
legend('F = 2, S = 3', 'F = 3, S = 3', 'F = 2, S = 6', 'F = 3, S = 6', 'Location', 'NorthEast')
xlabel('N')
ylabel('Recovery performance (median)')
grid on
title('Single-graph, multiple filters')

figure
semilogy(NN, plot_mean_data, 's-')
legend('F = 2, S = 3', 'F = 3, S = 3', 'F = 2, S = 6', 'F = 3, S = 6', 'Location', 'NorthEast')
xlabel('N')
ylabel('Recovery performance (mean)')
grid on
title('Single-graph, multiple filters')
