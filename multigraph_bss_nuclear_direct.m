function Z_hat = multigraph_bss_nuclear_direct(y, A, V, verbose)

if ~exist('verbose', 'var')
  verbose = false;
end

assert(length(A) == size(V, 3))
numGraphs = size(V, 3);

for i = 1:numGraphs-1
  assert(all(size(A{i}) == size(A{i+1})))
end

assert(size(A{1}, 1) == length(y))
N = size(A{1}, 1);
% It is assumed that the filters have the same order.
L = size(A{1}, 2)/N;

cvx_begin quiet
  variable Z(N, L, numGraphs);

  pho = 1;
  tau = 0.5;
  objective = 0;
  for i = 1:numGraphs
    objective = objective + pho*norm_nuc(Z(:, :, i)) + ...
      tau*sum(norms(Z(:, :, i)*diag(2.^(0:L-1)), 2, 2)); %#ok<*NODEF,*IDISVAR>
  end

  minimize(objective);

  subject to
    eq_constraint = 0;
    for i = 1:numGraphs
      eq_constraint = eq_constraint + V(:, :, i)*A{i}*vec(Z(:, :, i));
    end

    y == eq_constraint; %#ok<EQEFF>
cvx_end

if verbose
  svd_fmt_str = prepare_svd_fmt_str(size(Z(:, :, i), 2));

  a = 0;
  b = 0;
  for i = 1:numGraphs
    a = a + pho*norm_nuc(Z(:, :, i));
    b = b + tau*sum(norms(Z(:, :, i), 2, 2));
  end
  fprintf('%d %d\n', a, b)

  for i = 1:numGraphs
    fprintf(sprintf('svd(Z%d)=(%s)\n', i, svd_fmt_str), svd(Z(:, :, i)));
  end
end

Z_hat = Z;

end
