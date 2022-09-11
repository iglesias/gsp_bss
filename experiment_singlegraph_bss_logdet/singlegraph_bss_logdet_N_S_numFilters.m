function [success, recovery_performance, iters_to_solve] = singlegraph_bss_logdet_N_S_numFilters

num_simulations = 10;
verbose_bss_logdet = false;

params.L = 3;

NN = [60 100 140 180];
SS = [6 9];
NUM_FILTERS = [2 3];

for N = NN, for S = SS, for numFilters = NUM_FILTERS
  success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);
  recovery_performance = zeros(num_simulations, 1);

  params.N = N;
  params.S = S;
  params.numFilters = numFilters;

  parfor n = 1:num_simulations
    [truth, model, y] = singlegraph_svd_bss_gen_problem(params);
    [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);

    [UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
    Z_hat = zeros([size(Zsum_hat) numFilters]);
    for i = 1:numFilters
      Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
    end

    recovery_performance(n) = recovery_assessment_perms(truth.Z, Z_hat);
    if recovery_performance(n) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end
  end

  success_percent = sum(success)/num_simulations;
  fprintf('N%3d S%d L%d numFilters%d: success=%.2f\n', ...
          params.N, params.S, params.L, params.numFilters, ...
          success_percent);

%   save(sprintf('singlegraph_bss_logdet_N%d_S%d_L%d_numFilters%d', ...
%                params.N, params.S, params.L, params.numFilters));
end, end, end

end
