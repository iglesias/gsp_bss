function [truth, model, y] = bss_gen_problem_from_graphs(model, params)

if ~exist('params', 'var')
  params = struct;
end

numGraphs = length(model.G);

% Order of the filters (number of filter coefficients).
if isfield(params, 'L')
  L = params.L;
else
  L = 4;
end

% Number of nodes.
N = size(model.G(1).W, 1);
for i = 1:numGraphs
  [m, n] = size(model.G(i).W);
  assert(m == n)
  assert(N == size(model.G(i).W, 1))
end

% Number of non-zero input nodes.
if isfield(params, 'S')
  S = params.S;
else
  S = 1;
end

data_distribution = DataDistribution.Uniform;
shift_operator = ShiftOperator.Adjacency;

for i = 1:numGraphs
  % Graph Laplacian.
  model.G(i).L = diag(sum(model.G(i).W)) - model.G(i).W;

  switch shift_operator
  case ShiftOperator.Adjacency
    model.G(i).S = model.G(i).W;
  case ShiftOperator.Laplacian
    model.G(i).S = model.G(i).L;
  end

  [model.G(i).V, Lambda] = eig(model.G(i).S);
  model.G(i).U = inv(model.G(i).V);
  model.G(i).lambda = diag(Lambda);
  % disp(minmax(model.G(i).lambda'))
end

model.V = reshape([model.G.V], [N, N, numGraphs]);

% Filter coefficients.
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
  % disp(minmax(vec(model.Psi{i})'))
end

% Build filter matrices.
H = zeros(N, N * numGraphs);
for i = 1:numGraphs
  Hi = truth.h(1, i)*eye(N);
  for l = 1:L-1
    Hi = Hi + truth.h(l+1, i)*model.G(i).S^l;
  end

  H(:, N*(i-1)+1:N*i) = Hi;
end

% disp(minmax(H(:)'))

% Input.

truth.xSupport = zeros(numGraphs, S);
for i = 1:numGraphs
  truth.xSupport(i, :) = randperm(N, S);
end

truth.x = zeros(N, numGraphs);

switch data_distribution
case DataDistribution.Normal
  for i = 1:numGraphs
    truth.x(truth.xSupport(i, :), i) = randn(S, 1);
  end

case DataDistribution.Uniform
  for i = 1:numGraphs
    truth.x(truth.xSupport(i, :), i) = rand(S, 1);
  end
end

% Normalize input signals.
for i = 1:numGraphs
  truth.x(:, i) = truth.x(:, i) ./ norm(truth.x(:, i));
end

y = H*truth.x(:);

for i = 1:numGraphs
  model.A{i} = kr(model.Psi{i}', model.G(i).U')';
  % disp(minmax(vec(model.A{i})'))
end

for i = 1:numGraphs
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
end

end