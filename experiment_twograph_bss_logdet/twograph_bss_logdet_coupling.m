function twograph_bss_logdet_coupling

% Dependency to show progress of parfor:
% https://github.com/DylanMuir/ParforProgMon
% addpath ~/workspace/matlab/DylanMuir-ParforProgMon-9a1c257/

num_simulations = 1000;
verbose_multigraph_bss_logdet = false;

params.numGraphs = 2;

COUPLING = [0.0 0.4 0.7 0.9 0.95 1.0];
NN = [50 100];
LL = [3];
SS = [1 3];

for coupling = COUPLING, for N = NN, for L = LL, for S = SS
  tic
  success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);
  recovery_performance = zeros(num_simulations, 1);

  params.coupling = coupling;
  params.N = N;
  params.L = L;
  params.S = S;

  ppm = ParforProgMon(sprintf('coupling%03d ', coupling*100), num_simulations);
  parfor n = 1:num_simulations
    [truth, model, y] = multigraph_bss_gen_problem(params);
    [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, ...
                                          verbose_multigraph_bss_logdet);

    recovery_performance(n) = recovery_assessment(truth.Z, Z_hat);
    if recovery_performance(n) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end
    ppm.increment();
  end

  time = toc;
  success_percent = sum(success)/num_simulations;
  fprintf('coupling%03d N%d L%d S%d: success=%d time=%d\n', ...
          coupling*100, N, L, S, success_percent, time)

  save(sprintf('play_twograph_bss_logdet_coupling%03d_N%d_L%d_S%d', ...
               coupling*100, N, L, S));
end, end, end, end

end
