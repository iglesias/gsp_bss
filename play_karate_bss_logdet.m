function [err, err_sum, iters_to_solve] = play_karate_bss_logdet

% [err, err_sum, iters] = play_karate_bss_logdet;
% figure, hold on, plot(err), plot(err_sum), hold off

num_simulations = 100;
iters_to_solve = inf(num_simulations, 1);
params.S = 1;
err = zeros(num_simulations, 1);
err_sum = zeros(num_simulations, 1);
verbose = false;

parfor n = 1:num_simulations
  [truth, model, y] = karate_bss_gen_problem(params);
  [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose);
  
  
  [UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
  numFilters = length(truth.Z);
  Z_hat = zeros([size(Zsum_hat) numFilters]);
  for i = 1:numFilters
    Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
  end
  
  err(n) = recovery_assessment_perms(truth.Z, Z_hat);
  err_sum(n) = norm(truth.Zsum - sum(Z_hat, 3), 'fro') / norm(truth.Zsum, 'fro');
  iters_to_solve(n) = iter;
end

end