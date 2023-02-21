function alpha_1_directed_cycle

NN = 3:3:50;
alphas = work(NN);
plot_alphas(NN, alphas);
title('Vandermonde normalized')

NN = round(linspace(10, 1500, 19));
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

  S = 1; % [1, 2, 3] % With 3 the long plot doesn't end... within 13h
  assert(S <= N)

  L = 3;
  assert(L <= N)

  alphas(n) = compute_alpha_1(model, S, L, normalize_Psi);
end

end
