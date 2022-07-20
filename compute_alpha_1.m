function alpha_1 = compute_alpha_1(model, S, L, normalize_Psi)

numGraphs = length(model.G);
alpha_1 = zeros(numGraphs, 1);

for i = 1:numGraphs
  [model.G(i).V, Lambda] = eig(model.G(i).W);
  model.G(i).U = inv(model.G(i).V);

  N = size(model.G(i).W, 1);

  rho_U__S = rhof(model.G(i).U, S);

  model.G(i).lambda = diag(Lambda);
  model.Psi{i} = repmat(model.G(i).lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
  assert(N == size(model.Psi{i}, 1))

  if normalize_Psi
    [Psi_svd_L, ~, ~] = svd(model.Psi{i}, 0);
    model.Psi{i} = Psi_svd_L;
  else
    model.Psi{i} = model.Psi;
  end

  rho_Psi__L = rhof(model.Psi{i}, L);

  alpha_1(i) = 3/128 * 1/(rho_Psi__L * rho_U__S * log10(2*N*L*S));
end

end