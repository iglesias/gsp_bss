function brain_bss_logdet_S_numGraphs

% Dependency to show progress of parfor:
% https://github.com/DylanMuir/ParforProgMon
% addpath ~/workspace/matlab/DylanMuir-ParforProgMon-9a1c257/

num_simulations = 1000;
verbose_multigraph_bss_logdet = false;

params.L = 3;

NUM_GRAPHS = [2 3 4 5 6];
SS = [1 2 3 4 5];

for S = SS
  for numGraphs = NUM_GRAPHS
    tic
    success = zeros(num_simulations, 1);
    iters_to_solve = inf(num_simulations, 1);
    recovery_performance = zeros(num_simulations, 1);

    params.S = S;
    params.numGraphs = numGraphs;

%    ppm = ParforProgMon('Window', num_simulations);
    parfor n = 1:num_simulations
      [truth, model, y] = brain_bss_gen_problem(params);
      [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, ...
                                            verbose_multigraph_bss_logdet);

      recovery_performance(n) = recovery_assessment(truth.Z, Z_hat);
      if recovery_performance(n) < 1e-3
        success(n) = 1;
        iters_to_solve(n) = iter;
      end
%      ppm.increment();
    end

    time = toc;
    success_percent = sum(success)/num_simulations;
    fprintf('S%d numGraphs%d: success=%d time=%d\n', S, numGraphs, ...
            success_percent, time)

    save(sprintf('play_brain_bss_logdet_S%d_L%d_numGraphs%d', params.S, ...
         params.L, params.numGraphs));
  end
end

end
