NN = [60 100 140 180];
SS = zeros(1, 2);
NUM_FILTERS = zeros(1, 2);

% Data generated from 100 simulations with
% experiment_singlegraph_bss_logdet/singlegraph_bss_logdet_N_S_numFilters
% (it takes about 1.5h in Intel(R) Core™ i7-6700HQ CPU @ 2.60GHz and 16 GB
% of RAM)
%
% N 60 S6  numFilters2: success=0.44
% N 60 S6  numFilters3: success=0.01
% N 60 S9  numFilters2: success=0.03
% N 60 S9  numFilters3: success=0.00
% N100 S6  numFilters2: success=0.87
% N100 S6  numFilters3: success=0.44
% N100 S9  numFilters2: success=0.74
% N100 S9  numFilters3: success=0.03
% N140 S6  numFilters2: success=0.99
% N140 S6  numFilters3: success=0.94
% N140 S9  numFilters2: success=0.91
% N140 S9  numFilters3: success=0.53
% N180 S6  numFilters2: success=1.00
% N180 S6  numFilters3: success=0.96
% N180 S9  numFilters2: success=0.99
% N180 S9  numFilters3: success=0.78


% From the lines above
success_data = zeros(length(NN), length(SS), length(NUM_FILTERS));
% N=60 S=6
success_data(1,1,:) = [0.44 0.01];
% N=60 S=9
success_data(1,2,:) = [0.03 0.00];
% N=100 S=6
success_data(2,1,:) = [0.87 0.44];
% N=100 S=9
success_data(2,2,:) = [0.74 0.03];
% N=140 S=6
success_data(3,1,:) = [0.99 0.94];
% N=140 S=9
success_data(3,2,:) = [0.91 0.53];
% N=180 S=6
success_data(4,1,:) = [1.00 0.96];
% N=180 S=9
success_data(4,2,:) = [0.99 0.78];


num_plots = length(SS)*length(NUM_FILTERS);
plot_data = zeros(length(NN), num_plots);

for n = 1:length(NN)
  plot_idx = 1;
  for r = 1:2
    for s = 1:2
      plot_data(n, plot_idx) = success_data(n, s, r);
      plot_idx = plot_idx + 1;
    end
  end
end

figure
plot(NN, plot_data, 's--', 'LineWidth', 1.5)
legend('R = 2, S = 6', 'R = 2, S = 9', 'R = 3, S = 6', 'R = 3, S = 9', 'Location', 'SouthEast')
xlabel('N')
ylabel('Rate of successful recovery')
grid on