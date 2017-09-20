function [truth, model, y] = multigraph_bss_gen_problem

numGraphs = 3;
% Number of nodes.
N = 50;
% Edge existence probability.
p = 0.1;
% Adjacency matrices.
for i = 1:numGraphs
  model.G(i).W = generate_connected_ER(N, p);
  model.G(i).L = diag(sum(model.G(i).W))-model.G(i).W;

  assert(issymmetric(model.G(i).W))
  assert(issymmetric(model.G(i).L))

  [model.G(i).V, model.G(i).D] = eig(model.G(i).L);
  model.G(i).U = inv(model.G(i).V);
  model.G(i).lambda = diag(model.G(i).D);
end

model.V = reshape([model.G.V], [N, N, numGraphs]);

numFilterCoeffs = 2;
truth.h = zeros(numFilterCoeffs, numGraphs);
data_distibution = DataDistribution.Uniform;

switch data_distibution
case DataDistribution.Normal
  truth.h = randn(numFilterCoeffs, numGraphs);
  truth.h = truth.h / repmat(norms(truth.h, 2, 1), numFilterCoeffs, 1);

case DataDistribution.Uniform
  truth.h = rand(numFilterCoeffs, numGraphs);
end

for i = 1:numGraphs
  Psi{i} = repmat(model.G(i).lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);
end

% Build filter matrices.
H = zeros(N, N * numGraphs);
for i = 1:numGraphs
  Hi = truth.h(1, i)*eye(N);
  for l = 1:numFilterCoeffs-1
    Hi = Hi + truth.h(l+1, i)*model.G(i).L^l;
  end

  H(:, N*(i-1)+1:N*i) = Hi;
end

% Input.
% Number of non-zero input nodes.
S = 1;

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

switch data_distibution
case DataDistribution.Normal
  for i = 1:numGraphs
    truth.x(truth.xSupport(i, :), i) = randn(S, 1);
    truth.x(:, i) = truth.x(:, i) / norm(truth.x(:, i));
  end

case DataDistribution.Uniform
  for i = 1:numGraphs
    truth.x(truth.xSupport(i, :), i) = rand(S, 1);
  end
end

y = H*truth.x(:);

for i = 1:numGraphs
  model.A{i} = kr(Psi{i}', model.G(i).U')';
end

for i = 1:numGraphs
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
end

end
