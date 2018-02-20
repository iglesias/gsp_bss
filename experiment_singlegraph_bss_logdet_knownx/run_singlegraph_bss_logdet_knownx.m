function run_singlegraph_bss_logdet_knownx

params.N = 50;
params.L = 3;
work(params)

end

function work(params)

num_simulations = 100;
verbose_optimizer = false;
KNOWN_RATIO = [0.5 0.625 0.75];
PP = [3 4 1];
SS = [1 2 3 4 5];

for known_ratio = KNOWN_RATIO, for P = PP, for S = SS
  tic
  success = zeros(num_simulations, 1);
  sum_success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);
  recovery_performance = zeros(num_simulations, 1);
  sum_recovery_performance = zeros(num_simulations, 1);
  params.numFilters = P;
  params.S = S;

  parfor n = 1:num_simulations
    [truth, model, y] = singlegraph_bss_gen_problem(params);
    [Z_hat, iter, nonzero_idxs, zero_idxs] = singlegraph_bss_logdet_knownx(truth.x, y, model.A, model.G.V, known_ratio, verbose_optimizer);

    parsave(sprintf('20180125/singlegraph_bss_logdet_knownx_%s', ...
                    randomstring(20)), ...
            truth, model, y, Z_hat, iter, nonzero_idxs, zero_idxs)

    recovery_performance(n) = recovery_assessment(truth.Z, Z_hat);
    if recovery_performance(n) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end

    sum_recovery_performance(n) = norm(truth.Zsum - sum(Z_hat, 3), 'fro') / ...
                                  norm(truth.Zsum, 'fro');
    if sum_recovery_performance(n) < 1e-3
      sum_success(n) = 1;
    end

  end

  time = toc;
  success_ratio = sum(success)/num_simulations;
  sum_success_ratio = sum(sum_success)/num_simulations;
  fprintf('known=%.2f N%3d S%d L%d numFilters%d: success=%.2f sum_success=%.2f time=%5d\n', ...
          known_ratio, params.N, params.S, params.L, params.numFilters, ...
          success_ratio, sum_success_ratio, time)

%  save(sprintf('run_singlegraph_bss_logdet_knownx_known%03.0f_N%03d_S%d_L%d_numFilters%d', ...
%               known_ratio*100, params.N, params.S, params.L, params.numFilters));
end, end, end

end
