% Dependency to show progress of parfor:
% https://github.com/DylanMuir/ParforProgMon
% addpath ~/workspace/matlab/DylanMuir-ParforProgMon-9a1c257/

num_simulations = 100;
verbose_multigraph_bss_logdet = false;
NUM_NODES = [25 50 100];
TIME = zeros(1, length(NUM_NODES));
SUCCESS = zeros(1, length(NUM_NODES));

for m = 1:length(NUM_NODES)
  tic
  success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);

%  ppm = ParforProgMon(sprintf('num_nodes=%d ', NUM_NODES(m)), num_simulations);
  for n = 1:num_simulations
    params.N = NUM_NODES(m);
    [truth, model, y] = multigraph_bss_gen_problem(params);
    [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, ...
                                  verbose_multigraph_bss_logdet);
    if recovery_assessment(truth.Z, Z_hat) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end
%    ppm.increment();
  end

  SUCCESS(m) = sum(success)/num_simulations;
  TIME(m) = toc;
end

display(NUM_NODES)
display(SUCCESS)
display(TIME)
