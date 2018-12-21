clearvars

% Use problem generator just to generate a graph filter.
params.L = 1;
params.N = 20;
params.numFilters = 5;
[truth, model, y] = singlegraph_bss_gen_problem(params);

% Arbritrary signal at the output of the filter.
N = size(y, 1);
x = randn(N, 1);

%%
S = 4;
NK = nchoosek(1:N, S);

min_err = inf;

err = zeros(size(NK, 1), 1);
for i = 1:size(NK, 1)
  Islice = zeros(N, S);
  Islice(NK(i,:), :) = eye(S);
  cvx_begin quiet
    variable stilde(S, 1);
    minimize( norm(x - truth.H*Islice*stilde, 2) );
  cvx_end

  err(i) = norm(x - truth.H*Islice*stilde, 2);
  if err(i) < min_err
    min_err = err(i);
    best_stilde = stilde;
    best_Islice = Islice;
  end
end

s = best_Islice*best_stilde;
