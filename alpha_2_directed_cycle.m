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

  L = 2;
  S = 1; % [1 2]
  assert(L <= N)
  assert(S <= N)
  
  alphas(n) = compute_alpha_2(model, S, L, normalize_Psi);
end

end