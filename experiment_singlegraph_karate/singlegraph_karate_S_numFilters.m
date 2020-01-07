function singlegraph_karate_S_numFilters

verbose_bss_logdet = false;
num_simulations = 100;

params.L = 3;

SS = [2 3];
NUM_FILTERS = [2 3];
NOISE = [1e-6 1e-5 1e-4 1e-3 1e-2];

for S = SS, for numFilters = NUM_FILTERS, for noise = NOISE
  tic
  success = zeros(num_simulations, 1);
  iters_to_solve = inf(num_simulations, 1);
  recovery_performance = inf(num_simulations, 1);

  params.S = S;
  params.numFilters = numFilters;
  params.noise = noise;

  parfor n = 1:num_simulations
    [truth, model, y] = karate_bss_gen_problem(params);
    [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);
    
    if any(vec(isnan(Zsum_hat)))
      continue
    end
    
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

  time = toc;
  success_percent = sum(success)/num_simulations;
  fprintf('S%d numFilters%d noise%d: success=%.2f time=%5d\n', ...
          params.S, params.numFilters, params.noise, ...
          success_percent, time)

  save(sprintf('singlegraph_karate_S%d_numFilters%dnoise%d', ...
                params.S, params.numFilters, params.noise));
end, end, end

end
