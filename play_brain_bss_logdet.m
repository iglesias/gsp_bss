clearvars

num_simulations = 1;
verbose = false;

predictor = zeros(6, 6);

m = 1;
for brain_a = 1:6
  for brain_b = brain_a+1:6
    tic
    success = zeros(num_simulations, 1);
    iters_to_solve = inf(num_simulations, 1);

    params.brain_idxs = [brain_a brain_b];

%    ppm = ParforProgMon('Work', num_simulations);
    for n = 1:num_simulations
      [truth, model, y] = brain_bss_gen_problem(params);
%       [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, verbose);

%       if recovery_assessment(truth.Z, Z_hat) < 1e-3
%         success(n) = 1;
%         iters_to_solve(n) = iter;
%       end
%      ppm.increment();
    end

    fprintf('%d %d\n', brain_a, brain_b)
    predictor(brain_a, brain_b) = max(abs(model.G(1).U(:))) * max(abs(model.G(2).U(:))) * max(norms(model.Psi{1}', 2)'.^2) * max(norms(model.Psi{2}', 2)'.^2);
%     SUCCESS(m) = sum(success)/num_simulations
%     TIME(m) = toc
    m = m+1;
  end
end

% display(SUCCESS)
% display(TIME)
