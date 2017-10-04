function play_bss_logdet_svd

num_simulations = 100;
verbose_bss_logdet = false;
NUM_NODES = [20 40 60 80 100];
TIME = zeros(1, length(NUM_NODES));
SUCCESS = zeros(1, length(NUM_NODES));

for m = 1:length(NUM_NODES)
  tic
  success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);
  recovery_performance = zeros(num_simulations, 1);

%  ppm = ParforProgMon(sprintf('num_nodes=%d ', NUM_NODES(m)), num_simulations);
  parfor n = 1:num_simulations
    [truth, model, y] = bss_gen_problem(NUM_NODES(m));
    [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);

    [UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
    numFilters = length(truth.Z);
    Z_hat = zeros([size(Zsum_hat) numFilters]);
    for i = 1:numFilters
      Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
    end

    recovery_performance(n) = recovery_assessment_perms(truth.Z, Z_hat);
    if recovery_performance(n) < 1e-3
      success(n) = 1;
      iters_to_solve(n) = iter;
    end
%    ppm.increment();
  end

  SUCCESS(m) = sum(success)/num_simulations;
  TIME(m) = toc;

  save(sprintf('play_bss_logdet_svd_num_nodes=%d_S6_L3_F2', NUM_NODES(m)))
  fprintf('num_nodes=%d done!\n', NUM_NODES(m))
end

display(NUM_NODES)
display(SUCCESS)
display(TIME)

end
