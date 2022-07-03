function alpha_1_directed_cycle

NN = 3:3:50;
alphas = work(NN);
plot_alphas(NN, alphas);
title('Vandermonde normalized')

NN = 10:20:500;
normalize_Psi = false;
alphas = work(NN, normalize_Psi);
plot_alphas(NN, alphas, 'ro-')
title('Vandermonde w/o normalization')

end

function plot_alphas(NN, alphas, plotstr)

if nargin < 3
  plotstr = 'o-';
end

figure
plot(NN, alphas, plotstr, 'LineWidth', 1.5)
grid on
xlabel('N')
ylabel('\alpha_1')

end

function alphas = work(NN, normalize_Psi)

if nargin < 2
  normalize_Psi = true;
end

assert(isvector(NN))
alphas = zeros(length(NN), 1);

for n = 1:length(NN) 
  N = NN(n);

  % Adjacency matrix for the directed cycle.
  model.G.W = circshift(eye(N), [N 1]);

  [model.G.V, Lambda] = eig(model.G.W);
  U = inv(model.G.V);

  S = 1; % [1, 2, 3] % With 3 the long plot doesn't end... within 13h
  assert(S <= N)
  rho_U__S = -Inf;
  for ell = 1:N
    rho_U__S = max(rho_U__S, max(norms(nchoosek(U(ell, :), S), 2, 2).^2));
  end

  model.G.lambda = diag(Lambda);
  L = 3;
  assert(L <= N)
  model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
  assert(N == size(model.Psi, 1))

  if normalize_Psi
    [Psi_svd_L, ~, ~] = svd(model.Psi, 0);
    Psi = Psi_svd_L;
  else
    Psi = model.Psi;
  end

  rho_Psi__L = -Inf;
  for ell = 1:N
    rho_Psi__L = max(rho_Psi__L, max(norms(nchoosek(Psi(ell, :), L), 2, 2).^2));
  end

  alpha_1 = 3/128 * 1/(rho_Psi__L * rho_U__S * log10(2*N*L*S));
  alphas(n) = alpha_1;
end

end