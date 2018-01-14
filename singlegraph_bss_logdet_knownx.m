function [Z_hat, iter] = singlegraph_bss_logdet_knownx(x, y, A, V, known_ratio, verbose)

if ~exist('known_ratio', 'var')
  known_ratio = 0.5;
else
  if ~isreal(known_ratio) || ~isscalar(known_ratio) || ...
      known_ratio < 0 || known_ratio > 1
    error('known_ratio must be a scalar between 0 and 1.')
  end
end

if ~exist('verbose', 'var')
  verbose = false;
end

%% Parameters.

epsilon_rank = 5e-2;
epsilon_normx = 1e-3;
max_iter = 30;

%% Extract problem information from the inputs.

assert(size(A, 1) == length(y))
N = size(A, 1);
% It is assumed that the filters have the same order.
L = size(A, 2)/N;
B = V*A;
numFilters = size(x, 2);

%% Set up known values of the inputs.

num_known = round(known_ratio*N);
[nonzero_idxs, zero_idxs] = separate_known_idxs(x, num_known);

%% Optimization.

% Initialization.
flag = 1;
iter = 1;
Theta_old = repmat(eye(N), [1, 1, numFilters]);
Kappa_old = repmat(eye(L), [1, 1, numFilters]);
Z_old = zeros(N, L, numFilters);

% Majorization-minimization.
while (flag == 1 && iter <= max_iter)
  if verbose
    fprintf('iteration %d\n', iter)
  end

  %TODO vectorize
  wx = zeros(N, numFilters);
  for i = 1:numFilters
    wx(:, i) = 1./(sqrt(sum(abs(Z_old(:, :, i)).^2, 2)) + epsilon_normx);
  end

  cvx_begin quiet
    variable Z(N, L, numFilters);
    variable Theta(N, N, numFilters) symmetric;
    variable Kappa(L, L, numFilters) symmetric;

    objective = 0;
    for i = 1:numFilters
      objective = objective + ...
                  trace((Theta_old(:, :, i) + epsilon_rank*eye(N))\Theta(:, :, i)) + ...
                  trace((Kappa_old(:, :, i) + epsilon_rank*eye(L))\Kappa(:, :, i));
      objective = objective + wx(:, i)'*norms(Z(:, :, i), 2, 2);
    end

    minimize(objective);

    subject to
      for i = 1:numFilters
        [Theta(:, :, i) Z(:, :, i); Z(:, :, i)' Kappa(:, :, i)] == semidefinite(N+L);
      end

      y == B*vec(sum(Z, 3));

      for i = 1:numFilters
        Z(zero_idxs{i}, :, i) == 0;

        for j = 1:length(nonzero_idxs{i})
          for k = j+1:length(nonzero_idxs{i})
            jidx = nonzero_idxs{i}(j);
            kidx = nonzero_idxs{i}(k);
            if x(jidx, i)*x(kidx, i) == 0
              error('Unexpected, check for bug .')
            end
            if ~isempty(find(zero_idxs{i} == jidx, 1)) || ...
               ~isempty(find(zero_idxs{i} == kidx, 1))
              error('Unexpected, check for bug.')
            end
            x(jidx, i)*Z(kidx, :, i) == x(kidx, i)*Z(jidx, :, i);
          end
        end
      end
  cvx_end

  if verbose
    a = 0;
    b = 0;
    for i = 1:numFilters
      a = a + trace((Theta_old(:, :, i) + epsilon_rank*eye(N))\Theta(:, :, i)) + ...
              trace((Kappa_old(:, :, i) + epsilon_rank*eye(L))\Kappa(:, :, i));
      b = b + wx(:, i)'*norms(Z(:, :, i), 2, 2);
    end

    fprintf('%d %d\n', a, b)
  end

  if isempty(strfind(cvx_status, 'Solved'))
    fname = sprintf('failed_problem_%s_v%s', mfilename, ...
                    datestr(now, 'ddmmyyyyHHMMSS'));
    warning(sprintf('cvx_status not Solved, saving %s.', fname))
    save(fname)
    Z_hat = nan(size(Z));
    return
  end

  difference = 0;
  for i = 1:numFilters
    difference = difference + norm(Z(:, :, i) - Z_old(:, :, i), 'fro')/ ...
                              norm(Z_old(:, :, i), 'fro');
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
      if verbose
        fprintf('Convergence reached, cvx_status: %s.\n', cvx_status)
      end
      flag = 0;
    else
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

function [nonzero_idxs, zero_idxs] = separate_known_idxs(x, num_known)

[N, P] = size(x);

nonzero_idxs = cell(P, 1);
zero_idxs = cell(P, 1);

if num_known == 0
  return
end

if num_known > N
  error('The number of known input values cannot be larger than the signal dimension.')
end

known_idxs = zeros(num_known, P);
for i = 1:P
  known_idxs(:, i) = sort(randperm(N, num_known));
  nonzero_idxs{i} = known_idxs(find(x(known_idxs(:, i), i)), i);
  zero_idxs{i} = setdiff(known_idxs(:, i), nonzero_idxs{i});
end

end
