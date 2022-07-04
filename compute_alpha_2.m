function alpha_2 = compute_alpha_2(model, S, L, normalize_Psi)

[model.G.V, Lambda] = eig(model.G.W);
U = inv(model.G.V);

model.G.lambda = diag(Lambda);
N = size(model.G.W, 1);
assert(N == size(model.G.W, 2))

model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
assert(N == size(model.Psi, 1))

if normalize_Psi
  [Psi_svd_L, ~, ~] = svd(model.Psi, 0);
  Psi = Psi_svd_L;
else
  Psi = model.Psi;
end

kappa_U = kappaf(U, U, S, S);
kappa_Psi = kappaf(Psi, Psi, L, L);

rho_Psi1 = rhof(Psi, L);
rho_Psi2 = rhof(Psi, L);

rho_max = max(rhof(U, S), rhof(U, S));

% Number of sources R.
R = 2;
alpha_2 = 9/32 * 1/(rho_Psi1*rho_Psi2*rho_max + 1/R*kappa_U*kappa_Psi) * 1/(R^2*log10(2*N));

end