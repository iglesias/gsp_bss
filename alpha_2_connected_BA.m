function alpha_2s = alpha_2_connected_BA

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

alpha_2s = zeros(num_mc, 1);
for i = 1:num_mc
  alpha_2s(i) = work(N, S, L);
end

% disp(alpha_2s)

fprintf('min=%.4f median=%.4f max=%.4f\n', min(alpha_2s), median(alpha_2s), max(alpha_2s))

end

function alpha_2 = work(N, S, L)

model.G(1).W = generate_connected_BA(N, N/10);
model.G(2).W = generate_connected_BA(N, N/10);
normalize_Psi = true;
alpha_2 = compute_alpha_2(model, S, L, normalize_Psi);

end
