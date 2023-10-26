function alpha_2 = compute_alpha_2(model, S, L, normalize_Psi)

numGraphs = 2;
assert(length(model.G) == numGraphs)

for i = 1:numGraphs
  [model.G(i).V, Lambda] = eig(model.G.W);
  model.G(i).U = inv(model.G(i).V);
  model.G(i).lambda = diag(Lambda);
end

N = size(model.G(1).W, 1);
for i = 2:numGraphs
  assert(N == size(model.G(i).W, 1))
end

for i = 1:numGraphs
  model.Psi{i} = repmat(model.G(i).lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>' 
  assert(N == size(model.Psi{i}, 1))
  
  if normalize_Psi
    [Psi_svd_L, ~, ~] = svd(model.Psi{i}, 0);
    model.Psi{i} = Psi_svd_L;
  end
end

assert(numGraphs == 2)

U_1 = model.G(1).U;
U_2 = model.G(2).U;

Psi_1 = model.Psi{1};
Psi_2 = model.Psi{2};

kappa_U = kappaf(U_1, U_2, S, S);
kappa_Psi = kappaf(Psi_1, Psi_2, L, L);

rho_Psi1 = rhof(Psi_1, L);
rho_Psi2 = rhof(Psi_2, L);

rho_max = max(rhof(U_1, S), rhof(U_2, S));

R = numGraphs;
alpha_2 = (9/32 * 1/(rho_Psi1*rho_Psi2*rho_max + 1/R*kappa_U*kappa_Psi)) * 1/(R^2*log(N));

end