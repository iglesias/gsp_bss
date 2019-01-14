function [erdosrenyi, cycle]  =cc

% product of the norms of psi_2,l psi_1,l

L = 3;

N = 50;
p = 0.1;
for i = 1:2
  erdosrenyi(i).A = generate_connected_ER(N, p);
  [erdosrenyi(i).V, erdosrenyi(i).Lambda] = eig(erdosrenyi(i).A);
%  erdosrenyi(i).U = inv(erdosrenyi(i).V);
  erdosrenyi(i).Psi = make_Psi(diag(erdosrenyi(i).Lambda), L);
  erdosrenyi(i).Psi_norms = reshape(norms(erdosrenyi(i).Psi'), N, 1);
end

shift = [-1 1];
for i = 1:2
  cycle(i).A = circshift(eye(N), shift(i));
  [cycle(i).V, cycle(i).Lambda] = eig(cycle(i).A);
  %cycle(i).U = inv(cycle(i).V);
  cycle(i).Psi = make_Psi(diag(cycle(i).Lambda), L);
  cycle(i).Psi_norms = reshape(norms(cycle(i).Psi'), N, 1);
end

x = erdosrenyi(1).Psi_norms .* erdosrenyi(2).Psi_norms;
y = cycle(1).Psi_norms .* cycle(2).Psi_norms;

figure(1)
hist(x, 50)
title(sprintf('ER\nmedian=%d', median(x)))
xlabel('||\Psi_{1,l}|| ||\Psi_{2,l}||')
ylabel('Histogram count')

figure(2)
hist(y, 10)
title(sprintf('Cycle\nmedian=%d', median(y)))
xlabel('||\Psi_{1,l}|| ||\Psi_{2,l}||')
ylabel('Histogram count')

end

function Psi = make_Psi(lambda, L)

assert(isvector(lambda))
assert(isscalar(L))

N = numel(lambda);

Psi = repmat(lambda, 1, L).^repmat([0:L-1], N, 1);

end
