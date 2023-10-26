function alpha_3s = alpha_3_connected_BA

% Number of graph nodes.
N = 100;
% Number of graph realizations.
num_mc = 1e3;
% Number of non-zero nodes in the graph signal (sparsity strength).
S = 1;
assert(S <= N)
% Number of filter taps.
L = 1;
assert(L <= N)
fprintf('N=%d S=%d L=%d\n', N, S, L)

alpha_3s = zeros(num_mc, 1);
for i = 1:num_mc
  model.G.W = generate_connected_BA(N, N/10);
  normalize_Psi = true;
  alpha_3s(i) = compute_alpha_3(model, S, L, normalize_Psi);
end

fprintf('min=%.4f median=%.4f max=%.4f\n', min(alpha_3s), median(alpha_3s), max(alpha_3s))

end