function alpha_2_directed_cycle

NN = 3:3:50;
alphas = work(NN);
plot_alphas(NN, alphas);
title('Vandermonde normalized')

NN = 3:3:50;
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
ylabel('\alpha_2')

end

function alphas = work(NN, normalize_Psi)
% 5x in N -> 100x in time
% N  | time
% ----------
% 2  | 1e-4
% 10 | 1e-2
% 50 | 1e+0
% 250| 1e+2

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

  L = 2;
  S = 1; % [1 2]

  model.G.lambda = diag(Lambda);
  assert(L <= N)
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
  alphas(n) = alpha_2;
end

end