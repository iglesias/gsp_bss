function run_singlegraph_bss_logdet_knownx

%params.N = 50;
%params.L = 3;
%params.numFilters = 2;
%params.S = 5;
%work(params)
%clear params

params.N = 100;
params.L = 5;
params.numFilters = 4;
params.S = 5;
work(params)

end

function work(params)

num_simulations = 50;
verbose_optimizer = false;
KNOWN_RATIO = flip([.2 .4 .6 .8 1]);

for known_ratio = KNOWN_RATIO
  tic
  success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);
  recovery_performance = zeros(num_simulations, 1);

  parfor n = 1:num_simulations
    [truth, model, y] = singlegraph_bss_gen_problem(params);
    [Z_hat, iter] = singlegraph_bss_logdet_knownx(truth.x, y, model.A, model.G.V, known_ratio, verbose_optimizer);

    iter = randi(30);
    recovery_performance(n) = recovery_assessment(truth.Z, Z_hat);
    if recovery_performance(n) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end
  end

  time = toc;
  success_ratio = sum(success)/num_simulations;
  fprintf('known=%.2f N%3d S%d L%d numFilters%d: success=%.2f time=%5d\n', ...
          known_ratio, params.N, params.S, params.L, params.numFilters, ...
          success_ratio, time)

  save(sprintf('run_singlegraph_bss_logdet_knownx_known%03.0f_N%03d_S%d_L%d_numFilters%d', ...
               known_ratio*100, params.N, params.S, params.L, params.numFilters));
end

end
