function p = probability_connected_BA

% Number of graph nodes.
N = 100;
% Number of graph realizations.
num_mc = 1e1;
% Number of non-zero nodes in the graph signal (sparsity strength).
S = 1;
assert(S <= N)
% Number of filter taps.
L = 1;
assert(L <= N)
fprintf('N=%d S=%d L=%d\n', N, S, L)

p = zeros(num_mc, 1);
for i = 1:num_mc
  model.G(1).W = generate_connected_BA(N, N/10);
  model.G(2).W = generate_connected_BA(N, N/10);
  normalize_Psi = true;
  alpha_1 = compute_alpha_1(model, S, L, normalize_Psi);
  alpha_3 = compute_alpha_3(model, S, L, normalize_Psi);
  alpha = min([alpha_1; alpha_3]);
  p(i) = 1 - N^(-alpha+1);
end

hist(p)

end
