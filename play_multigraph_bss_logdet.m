% Dependency to show progress of parfor:
% https://github.com/DylanMuir/ParforProgMon
% addpath ~/workspace/matlab/DylanMuir-ParforProgMon-9a1c257/

N = 100;
verbose_multigraph_bss_logdet = false;
NUM_NODES = [25 50 100];
TIME = zeros(1, length(NUM_NODES));
SUCCESS = zeros(1, length(NUM_NODES));

for m = 1:length(NUM_NODES)
  tic
  success = zeros(N, 1);
  iters_to_solve = inf(N, 1);

  ppm = ParforProgMon(sprintf('num_nodes=%d ', NUM_NODES(m)), N);
  parfor n = 1:N
    [truth, model, y] = multigraph_bss_gen_problem(NUM_NODES(m));
    [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, ...
                                  verbose_multigraph_bss_logdet);
    if recovery_assessment(truth.Z, Z_hat) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end
    ppm.increment();
  end

  SUCCESS(m) = sum(success)/N;
  TIME(m) = toc;
end

display(NUM_NODES)
display(SUCCESS)
display(TIME)
