N = 50;
verbose_bss_logdet = false;
NUM_NODES = [50];
TIME = zeros(1, length(NUM_NODES));
SUCCESS = zeros(1, length(NUM_NODES));

for m = 1:length(NUM_NODES)
  tic
  success = zeros(N, 1);
  iters_to_solve = inf(N, 1);

  ppm = ParforProgMon(sprintf('num_nodes=%d ', NUM_NODES(m)), N);
  parfor n = 1:N
    [truth, model, y] = bss_gen_problem(NUM_NODES(m));
    [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);

    [UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
    numFilters = length(truth.Z);
    Z_hat = zeros([size(Zsum_hat) numFilters]);
    for i = 1:numFilters
      Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
    end

    if recovery_assessment_perms(truth.Z, Z_hat) < 1e-3
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
