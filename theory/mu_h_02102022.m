clearvars

L = 5;
%%%% Erdos-Renyi:
% N = 50;
% % Edge existence probability.
% p = 0.1;
% % Adjacency matrix.
% model.G.W = generate_connected_ER(N, p);
%
%%% Directed cycle:
% N = 50;
% model.G.W = circshift(eye(N), [N 1]);
%
%%%% Brain:
params.numGraphs = 1;
[~, model, ~] = brain_bss_gen_problem(params);
N = size(model.G.W, 1);
%
normalize_Psi = true;
assert(L <= N)

[model.G.V, Lambda] = eig(model.G.W);
model.G.U = inv(model.G.V);
U = model.G.U;

model.G.lambda = diag(Lambda);
model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
assert(N == size(model.Psi, 1))

if normalize_Psi
  [Psi_svd_L, ~, ~] = svd(model.Psi, 0);
  Psi = Psi_svd_L;
else
  Psi = model.Psi; %#ok<UNRCH>
end

numFilters = 100;

h1 = randn(L, numFilters);
% Heat Kernel filter coefficients.
h2 = exp(-(0.3*rand(numFilters, 1)+0.7)*(L:-1:1))';

% Normalize filter taps.
h1 = h1 ./ repmat(norms(h1, 2, 1), L, 1);
h2 = h2 ./ repmat(norms(h2, 2, 1), L, 1);

mu_h_squared_1 = N*max(abs(conj(Psi)*h1).^2);
mu_h_squared_2 = N*max(abs(conj(Psi)*h2).^2);

b = ceil(max([mu_h_squared_1 mu_h_squared_2]));

subplot(211)
histogram(mu_h_squared_1, linspace(1, b, 50))
yl = ylim;
hold on
plot(mean(mu_h_squared_1)*[1 1], [0 yl(2)], 'lineWidth', 2)
hold off
title('{\mu_h}^2 with Gaussian filters')

subplot(212)
hist(mu_h_squared_2, linspace(1, b, 50))
yl = ylim;
hold on
plot(mean(mu_h_squared_2)*[1 1], [0 yl(2)], 'lineWidth', 2)
hold off
title('{\mu_h}^2 with heat kernel filters')
xlabel(abs(mean(mu_h_squared_2)-mean(mu_h_squared_1)))