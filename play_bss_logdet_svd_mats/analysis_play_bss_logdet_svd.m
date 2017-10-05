NN = [20 40 60 80 100];
SS = [3 6];
FF = [2 3];

num_plots = length(SS)*length(FF);
plot_data = zeros(length(NN), num_plots);
fname_pattern = 'play_bss_logdet_svd_num_nodes=%d_S%d_L3_F%d';

for n = 1:length(NN)
  plot_idx = 1;
  for S = SS
    for F = FF
      fname = sprintf(fname_pattern, NN(n), S, F);
      fprintf('Reading %s\n', fname);
      load(fname)

      plot_data(n, plot_idx) = median(recovery_performance);

      plot_idx = plot_idx + 1;
    end
  end
end

semilogy(NN, plot_data, 's-')
legend('F = 2, S = 3', 'F = 3, S = 3', 'F = 2, S = 6', 'F = 3, S = 6', 'Location', 'East')
xlabel('N')
ylabel('Recovery performance')
grid on
title('Single-graph, multiple filters')
