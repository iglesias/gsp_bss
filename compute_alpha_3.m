function alpha_3 = compute_alpha_3(model, S, L, normalize_Psi)

numGraphs = length(model.G);
alpha_3 = zeros(numGraphs, 1);

h = randn(L, 1);
% Normalize filter taps.
h = h / norm(h);

for i = 1:numGraphs
  [model.G(i).V, Lambda] = eig(model.G(i).W);
  model.G(i).U = inv(model.G(i).V);

  N = size(model.G(i).W, 1);

  model.G(i).lambda = diag(Lambda);
  model.Psi{i} = repmat(model.G(i).lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
  assert(N == size(model.Psi{i}, 1))

  if normalize_Psi
    [Psi_svd_L, ~, ~] = svd(model.Psi{i}, 0);
    model.Psi{i} = Psi_svd_L;
  else
    model.Psi{i} = model.Psi;
  end

  mu_max_squared = max(N/L*max(norms(model.Psi{i}, 2, 2).^2));
  %   display(mu_max_squared)
  mu_h_squared = N*max(max(abs(conj(model.Psi{i})*h).^2) / norm(h)^2);
  %   display(mu_h_squared)

  R = 2;

  % After review (decoupled model):
  rho_max = rhof(model.G(i).U, 1);
  alpha_3(i) = N / (R*(rho_max + 1)^2 * max(mu_max_squared*L, mu_h_squared*S^2) * (log(N))^2) - log(R);
end
  
end