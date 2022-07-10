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
  alpha_1s(i) = work(N, S, L, normalize_Psi);
end

% disp(alpha_1s)

fprintf('min=%.4f median=%.4f max=%.4f\n', min(alpha_1s), median(alpha_1s), max(alpha_1s))

end

function alpha_1 = work(N, S, L, normalize_Psi)

if nargin < 3
  normalize_Psi = true;
end

model.G.W = generate_connected_BA(N, N/10);

[model.G.V, Lambda] = eig(model.G.W);
U = inv(model.G.V);

rho_U__S = rhof(U, S);

model.G.lambda = diag(Lambda);
model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
assert(N == size(model.Psi, 1))

if normalize_Psi
  [Psi_svd_L, ~, ~] = svd(model.Psi, 0);
  Psi = Psi_svd_L;
else
  Psi = model.Psi;
end

rho_Psi__L = rhof(Psi, L);

alpha_1 = 3/128 * 1/(rho_Psi__L * rho_U__S * log10(2*N*L*S));

end
