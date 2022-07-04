function alpha_2s = alpha_2_connected_BA

% Number of graph nodes.
N = 50;
% Number of graph realizations.
num_mc = 50;
% Number of non-zero nodes in the graph signal (sparsity strength).
S = 1;
assert(S <= N)
% Number of filter taps.
L = 1;
assert(L <= N)

alpha_2s = zeros(num_mc, 1);
for i = 1:num_mc
  alpha_2s(i) = work(N, S, L);
end

% disp(alpha_2s)

end

function alpha_2 = work(N, S, L)

% Adjacency matrix of a Barabábasi-Albert (BA) graph.
model.G.W = generate_connected_BA(N, N/10);
normalize_Psi = true;
alpha_2 = compute_alpha_2(model, S, L, normalize_Psi);

end