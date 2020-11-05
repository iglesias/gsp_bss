clearvars

num_simulations = 10;
verbose = false;

predictor = zeros(6, 6);
probs = zeros(6, 6);

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
      [Z_hat, iter] = multigraph_bss_logdet(y, model.A, model.V, verbose);

      if recovery_assessment(truth.Z, Z_hat) < 1e-3
        success(n) = 1;
        iters_to_solve(n) = iter;
      end
%      ppm.increment();
    end

    % rho(U, k=1) for U with size(U)=[N,N],       rho(U, k=1)   = max_{ele=1:N} max_{j=1:N} abs(A(ele,j))
    % rho(Psi, k=2) for Psi with size(Psi)=[N,2], rho(Psi, k=2) = max_{ele=1:N} norm(Psi(ele,:), 2)^2
    predictor(brain_a, brain_b) = max(abs(model.G(1).U(:))) * max(abs(model.G(2).U(:))) * max(norms(model.Psi{1}', 2)'.^2) * max(norms(model.Psi{2}', 2)'.^2);
    
    % kappa(U1, U2, k1=1, k2=1) for U1, U2 with size(U1)=size(U2)=[N,N]
    % kappa(U1, U2, k1=1, k2=1) = max_{ele=1:N} max_{i=1:N} max_{j=1:N} norm( U1(ele,i) * U2(ele,j)' )
    % kappa(Psi1, Psi2, k1=2, k2=2) for Psi1, Psi2 with size(Psi1)=size(Psi2)=[N,2]
    % kappa(Psi1, Psi2, k1=2, k2=2) = max_{ele=1:N} norm( Psi1(ele,:) * Psi2(ele,:)' )
    
    probs(brain_a, brain_b) = sum(success)/num_simulations;

   fprintf('%d %d\n', brain_a, brain_b)
    SUCCESS(m) = sum(success)/num_simulations
    TIME(m) = toc
    m = m+1;
  end
end

% display(SUCCESS)
% display(TIME)

% probs = get_brain_pairs_success;
analyse_brain_pairs(predictor, probs);
