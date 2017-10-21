function  [truth, model, y] = brain_bss_gen_problem(params)

load('data/brain_data_66', 'CC')
assert(size(CC, 1) == size(CC, 2))
% Number of nodes.
N = size(CC, 2);

if ~exist('params', 'var')
  params = struct;
end

if isfield(params, 'brain_idxs')
  brain_idxs = params.brain_idxs;
else
  if isfield(params, 'numGraphs')
    if params.numGraphs > size(CC, 3)
      error('Maximum number of graphs is %d (%d provided)', size(CC, 3), ...
            params.numGraphs);
    end

    brain_idxs = randperm(size(CC, 3), params.numGraphs);
  else
    brain_idxs = randperm(size(CC, 3), 2);
  end
end

numGraphs = length(brain_idxs);

% Order of the filters (number of filter coefficients).
if isfield(params, 'L')
  L = params.L;
else
  L = 2;
end

% Number of non-zero input nodes.
if isfield(params, 'S')
  S = params.S;
else
  S = 1;
end

data_distribution = DataDistribution.Normal;
shift_operator = ShiftOperator.Adjacency;

for i = 1:numGraphs
  % Adjacency matrix.
  brain_graph = reshape(CC(:,:,brain_idxs(i)), N, N)*100;
  model.G(i).W = double(brain_graph >= min(max(brain_graph)));
  % Graph Laplacian.
  model.G(i).L = diag(sum(model.G(i).W))-model.G(i).W;

  assert(issymmetric(model.G(i).W))
  assert(issymmetric(model.G(i).L))

  switch shift_operator
  case ShiftOperator.Adjacency
    model.G(i).S = model.G(i).W;
  case ShiftOperator.Laplacian
    model.G(i).S = model.G(i).L;
  end

  [model.G(i).V, Lambda, model.G(i).U] = eig(model.G(i).S);
  model.G(i).lambda = diag(Lambda);
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
  model.A{i} = kr(model.Psi{i}', model.G(i).U)';
end

for i = 1:numGraphs
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
end

end
