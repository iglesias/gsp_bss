function [Z_hat, iter] = bss_logdet_jointsum(y, A, V, verbose)

if ~exist('verbose', 'var')
  verbose = false;
end

%% Parameter definition.

epsilon_rank = 5e-2;
epsilon_normx = 1e-3;
max_iter = 50;

%% Blind source separation using rank mininimization: log-det surrogate

assert(size(A, 1) == length(y))
N = size(A, 1);
% It is assumed that the filters have the same order.
L = size(A, 2)/N;
B = V*A;

% Initialization
flag = 1;
iter = 1;
Theta_old = eye(N);
Kappa_old = eye(L);

% LS solution for initializating w.
Z_old = reshape(pinv(B)*y, N, L);

% Majorization-minimization
while (flag == 1 && iter <= max_iter)
  if verbose
    fprintf('iteration %d\n', iter)
  end

  wx = 1./(sqrt(sum(abs(Z_old).^2, 2)) + epsilon_normx);

  cvx_begin quiet
    variable Z(N, L);
    variable Theta(N, N) symmetric;
    variable Kappa(L, L) symmetric;

    pho = 1;
    tau = 0.25;
    minimize(pho*(trace((Theta_old + epsilon_rank*eye(N))\Theta) + ...
                  trace((Kappa_old + epsilon_rank*eye(L))\Kappa)) + ...
             tau*wx'*norms(Z, 2, 2));

    subject to
      [Theta Z; Z' Kappa] == semidefinite(N+L);
      y == B*Z(:);
  cvx_end

  if verbose
    a = pho*(trace((Theta_old + epsilon_rank*eye(N))\Theta) + ...
             trace((Kappa_old + epsilon_rank*eye(L))\Kappa));
    b = tau*wx'*norms(Z, 2, 2);
    fprintf('%d %d\n', a, b)
  end

  if isempty(strfind(cvx_status, 'Solved'))
    fname = sprintf('failed_problem_bss_logdet_jointsum_v%s', ...
                    datestr(now, 'ddmmyyyyHHMMSS'));
    warning(sprintf('cvx_status not Solved, saving %s.', fname))
    save(fname)
    Z_hat = nan(size(Z));
    return
  end

  difference = norm(Z - Z_old, 'fro')/norm(Z_old, 'fro');

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
      Theta_old = Theta;
      Kappa_old = Kappa;
      Z_old = Z;
      iter = iter + 1;
    end
  end
end

Z_hat = Z;

end
