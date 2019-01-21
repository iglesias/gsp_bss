function [success, iters_to_solve] = assess_recovery_mixed_graphs

verbose_multigraph_bss_logdet = false;

num_simulations = 100;
success = zeros(num_simulations, 1);
iters_to_solve = inf(num_simulations, 1);

N = 50;
p = 0.8;

for n = 1:num_simulations
  model.G(1).W = generate_connected_ER(N, p);
  model.G(2).W = circshift(eye(N), 1);

  [truth, model, y] = bss_gen_problem_from_graphs(model);
  [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, verbose_multigraph_bss_logdet);

  if recovery_assessment(truth.Z, Z_hat) < 1e-3
    success(n) = 1;
    iters_to_solve(n) = iter;
  end
end

fprintf('%d\n', sum(success)/num_simulations);

end
