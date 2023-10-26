function alpha_1s = alpha_1_connected_BA

% Number of graph nodes.
N = 50;
% Number of graph realizations.
num_mc = 1e1;
% Number of non-zero nodes in the graph signal (sparsity strength).
S = 1;
assert(S <= N)
% Number of filter taps.
L = 1;
assert(L <= N)
fprintf('N=%d S=%d L=%d\n', N, S, L)

alpha_1s = zeros(num_mc, 1);
for i = 1:num_mc
  normalize_Psi = true;
  model.G.W = generate_connected_BA(N, N/10);
  alpha_1s(i) = compute_alpha_1(model, S, L, normalize_Psi);
end

% disp(alpha_1s)

fprintf('min=%.4f median=%.4f max=%.4f\n', min(alpha_1s), median(alpha_1s), max(alpha_1s))

end