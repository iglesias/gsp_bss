function [Z_hat, iter] = multigraph_bss_logdet(y, A, V, verbose)

if ~exist('verbose', 'var')
  verbose = false;
end

%% Parameter definition.

epsilon_rank = 5e-2;
epsilon_normx = 1e-3;
max_iter = 5;

assert(length(A) == size(V, 3))
numGraphs = size(V, 3);

%% Blind source separation using rank mininimization: log-det surrogate
for i = 1:numGraphs-1
  assert(all(size(A{i}) == size(A{i+1})))
end

assert(size(A{1}, 1) == length(y))
N = size(A{1}, 1);
% It is assumed that the filters have the same order.
L = size(A{1}, 2)/N;

% Initialization
flag = 1;
iter = 1;
Theta_old = repmat(eye(N), [1, 1, numGraphs]);
Kappa_old = repmat(eye(L), [1, 1, numGraphs]);

% Find initial values.
Z_old = multigraph_bss_nuclear_direct(y, A, V, verbose);

% Majorization-minimization
while (flag == 1 && iter <= max_iter)
  if verbose
    fprintf('iteration %d\n', iter)
  end

  %TODO vectorize
  wx = zeros(N, numGraphs);
  for i = 1:numGraphs
    wx(:, i) = 1./(sqrt(sum(abs(Z_old(:, :, i)).^2, 2)) + epsilon_normx);
  end

  cvx_begin quiet
    variable Z(N, L, numGraphs);
    variable Theta(N, N, numGraphs) symmetric;
    variable Kappa(L, L, numGraphs) symmetric;

    pho = 1;
    tau = 0.5;
    objective = 0;

    for i = 1:numGraphs
      objective = objective + ...
        pho*(trace((Theta_old(:, :, i) + epsilon_rank*eye(N))\Theta(:, :, i)) + ...
        trace((Kappa_old(:, :, i) + epsilon_rank*eye(L))\Kappa(:, :, i))) + ...
        tau*wx(:, i)'*norms(Z(:, :, i), 2, 2);
    end

    minimize(objective);

    subject to
      for i = 1:numGraphs
        [Theta(:, :, i) Z(:, :, i); Z(:, :, i)' Kappa(:, :, i)] == semidefinite(N+L);
      end

      eq_constraint = 0;
      for i = 1:numGraphs
        eq_constraint = eq_constraint + V(:, :, i)*A{i}*vec(Z(:, :, i));
      end

      y == eq_constraint;
  cvx_end

  if verbose
    a = 0;
    b = 0;
    for i = 1:numGraphs
      a = a + pho*(trace((Theta_old(:, :, i) + epsilon_rank*eye(N))\Theta(:, :, i)) + ...
        trace((Kappa_old(:, :, i) + epsilon_rank*eye(L))\Kappa(:, :, i)));
      b = b + tau*wx(:, i)'*norms(Z(:, :, i), 2, 2);
    end

    fprintf('%d %d\n', a, b)
  end

  if isempty(strfind(cvx_status, 'Solved'))
    fname = sprintf('failed_problem_multigraph_bss_logdet_v%s', ...
                    datestr(now, 'ddmmyyyyHHMMSS'));
    warning(sprintf('cvx_status not Solved, saving %s.', fname))
    save(fname)
    Z_hat = nan(size(Z));
    return
  end

  difference = 0;
  for i = 1:numGraphs
    difference = difference + norm(Z(:, :, i) - Z_old(:, :, i), 'fro')/norm(Z(:, :, i), 'fro');
  end

  if ~isempty(strfind(cvx_status, 'Infeasible'))
    % Stop the algorithm.
    if verbose
      fprintf('Infeasible cvx_status.\n')
    end
    Z_hat = nan(size(Z));
    return
  else
    if difference < 1e-4
      % Converged.
      if verbose
        fprintf('Convergence reached, cvx_status: %s.\n', cvx_status)
      end
      flag = 0;
    else
      % Did not converge.
      if verbose
        fprintf('Convergence NOT reached, difference=%d.\n', difference)
      end
      Z_old = Z;
      Theta_old = Theta;
      Kappa_old = Kappa;
      iter = iter + 1;
    end
  end
end

Z_hat = Z;

end
