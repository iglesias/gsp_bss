function [truth, model, y] = multigraph_bss_gen_problem(num_nodes)

numGraphs = 2;
% Number of filter coefficients.
L = 2;
data_distribution = DataDistribution.Uniform;
% Number of non-zero input nodes.
S = 1;

% Number of nodes.
if exist('num_nodes', 'var')
  N = num_nodes;
else
  N = 50;
end

% Edge existence probability.
p = 0.1;
% Adjacency matrices.
for i = 1:numGraphs
  model.G(i).W = generate_connected_ER(N, p);
  model.G(i).L = diag(sum(model.G(i).W))-model.G(i).W;

  assert(issymmetric(model.G(i).W))
  assert(issymmetric(model.G(i).L))

  [model.G(i).V, model.G(i).D] = eig(model.G(i).W);
  model.G(i).U = inv(model.G(i).V);
  model.G(i).lambda = diag(model.G(i).D);
end

model.V = reshape([model.G.V], [N, N, numGraphs]);

truth.h = zeros(L, numGraphs);

switch data_distribution
case DataDistribution.Normal
  truth.h = randn(L, numGraphs);
  truth.h = truth.h ./ repmat(norms(truth.h, 2, 1), L, 1);

case DataDistribution.Uniform
  truth.h = rand(L, numGraphs);
  truth.h = truth.h ./ repmat(norms(truth.h, 2, 1), L, 1);
end

for i = 1:numGraphs
  model.Psi{i} = repmat(model.G(i).lambda, 1, L).^repmat([0:L-1], N, 1);
end

% Build filter matrices.
H = zeros(N, N * numGraphs);
for i = 1:numGraphs
  Hi = truth.h(1, i)*eye(N);
  for l = 1:L-1
    Hi = Hi + truth.h(l+1, i)*model.G(i).W^l;
  end

  H(:, N*(i-1)+1:N*i) = Hi;
end

% Input.

truth.xSupport = zeros(numGraphs, S);
for i = 1:numGraphs
  truth.xSupport(i, :) = randperm(N, S);
end

%TODO generalize intersection from 2 graphs to numGraphs.
% if isempty(intersect(truth.x1Support, truth.x2Support))
%   fprintf('Intersection of the supports of the inputs empty.\n')
% else
%   fprintf('Intersection of the supports of the inputs non-empty.\n')
% end

truth.x = zeros(N, numGraphs);

switch data_distribution
case DataDistribution.Normal
  for i = 1:numGraphs
    truth.x(truth.xSupport(i, :), i) = randn(S, 1);
    truth.x(:, i) = truth.x(:, i) ./ norm(truth.x(:, i));
  end

case DataDistribution.Uniform
  for i = 1:numGraphs
    truth.x(truth.xSupport(i, :), i) = rand(S, 1);
    truth.x(:, i) = truth.x(:, i) ./ norm(truth.x(:, i));
  end
end

y = H*truth.x(:);

for i = 1:numGraphs
  model.A{i} = kr(model.Psi{i}', model.G(i).U')';
end

for i = 1:numGraphs
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
end

end
