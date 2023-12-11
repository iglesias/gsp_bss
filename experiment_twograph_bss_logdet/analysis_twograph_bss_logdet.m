COUPLING = [0.0 0.4 0.7 0.9 0.95 1.0];
NN = [100 50];
LL = [3];
SS = [1 3];

num_plots = length(NN)*length(LL)*length(SS);
plot_mean_data = zeros(length(COUPLING), num_plots);
plot_median_data = zeros(length(COUPLING), num_plots); 
plot_success_data = zeros(length(COUPLING), num_plots);
fname_pattern = 'play_twograph_bss_logdet_coupling%03d_N%d_L%d_S%d';
legend_str = {};

for c = 1:length(COUPLING)
  plot_idx = 1;
  for S = SS, for N = NN, for L = LL,
    fname = sprintf(fname_pattern, COUPLING(c)*100, N, L, S);
    fprintf('Reading %s\n', fname);
    load(fname, 'recovery_performance', 'success_percent')

%     if N == 50 && L == 3 && S == 3
%       figure
%       hist(recovery_performance)
%       title(sprintf('Coupling = %03d, N = %d, L = %d, S = %d', ...
%                     COUPLING(c)*100, N, L, S));
%       xlabel('Recovery performance')
%     end

%     plot_mean_data(c, plot_idx) = mean(recovery_performance);
%     plot_median_data(c, plot_idx) = median(recovery_performance);
    plot_success_data(c, plot_idx) = sum((recovery_performance/2)<1e-3)/length(recovery_performance);
    
    legend_str{plot_idx} = sprintf('N = %3d, S = %d', N, S);
    plot_idx = plot_idx + 1;
  end, end, end
end

% figure
% semilogy(COUPLING, plot_median_data, 's--')
% legend(legend_str)
% xlabel('Coupling', 'FontSize', 14)
% ylabel('RMSE (median)', 'FontSize', 14)
% grid on

figure
plot(COUPLING, plot_success_data, 's--', 'LineWidth', 1.5)
legend(legend_str, 'Location', 'SouthWest')
xlabel('Coupling', 'FontSize', 14)
ylabel('Rate of successful recovery', 'FontSize', 14)
set(gca, 'FontSize', 14)
grid on
