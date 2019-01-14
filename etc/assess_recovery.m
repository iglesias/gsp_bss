function [success, iters_to_solve] = assess_recovery

verbose_multigraph_bss_logdet = false;

num_simulations = 100;
success = zeros(num_simulations, 1);
iters_to_solve = inf(num_simulations, 1);

for n = 1:num_simulations
  [truth, model, y] = multigraph_bss_gen_problem;
  [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, verbose_multigraph_bss_logdet);

  if recovery_assessment(truth.Z, Z_hat) < 1e-3
    success(n) = 1;
    iters_to_solve(n) = iter;
  end
end

fprintf('%d\n', sum(success)/num_simulations);

end
