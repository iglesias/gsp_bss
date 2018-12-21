function [truth, model, y] = singlegraph_bss_gen_problem(params)

if ~exist('params', 'var')
  params = struct;
end

if ~isfield(params, 'verbose')
  params.verbose = false;
end

if isfield(params, 'numFilters')
  numFilters = params.numFilters;
else
  numFilters = 3;
end

% Order of the filters (number of filter coefficients).
if isfield(params, 'L')
  L = params.L;
else
  L = 2;
end

% Number of nodes.
if isfield(params, 'N')
  N = params.N;
else
  N = 50;
end

% Number of non-zero input nodes.
if isfield(params, 'S')
  S = params.S;
else
  S = 1;
end

data_distribution = DataDistribution.Normal;
shift_operator = ShiftOperator.Adjacency;

% Edge existence probability.
p = 0.1;
% Adjacency matrix.
model.G.W = generate_connected_ER(N, p);
% Graph Laplacian.
model.G.L = diag(sum(model.G.W))-model.G.W;

assert(issymmetric(model.G.W))
assert(issymmetric(model.G.L))

switch shift_operator
case ShiftOperator.Adjacency
  model.G.S = model.G.W;
case ShiftOperator.Laplacian
  model.G.S = model.G.L;
end

[model.G.V, Lambda, model.G.U] = eig(model.G.S);
model.G.lambda = diag(Lambda);

% Filter coefficients.
truth.h = zeros(L, numFilters);

switch data_distribution
case DataDistribution.Normal
  truth.h = randn(L, numFilters);
  truth.h = truth.h ./ repmat(norms(truth.h, 2, 1), L, 1);

case DataDistribution.Uniform
  truth.h = rand(L, numFilters);
  truth.h = truth.h ./ repmat(norms(truth.h, 2, 1), L, 1);
end

model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1);

% Build filter matrices.
truth.H = zeros(N, N * numFilters);
for i = 1:numFilters
  Hi = truth.h(1, i)*eye(N);
  for l = 1:L-1
    Hi = Hi + truth.h(l+1, i)*model.G.S^l;
  end

  Hi*randn(N, 1)
  truth.H(:, N*(i-1)+1:N*i) = Hi;
end

% Input signals.
truth.xSupport = zeros(numFilters, S);
for i = 1:numFilters
  truth.xSupport(i, :) = randperm(N, S);
end

if params.verbose
  found_overlap = false;
  i = 1;
  while i <= numFilters && ~found_overlap
    for j = i+1:numFilters
      if ~isempty(intersect(truth.xSupport(i, :), truth.xSupport(j, :)))
        found_overlap = true;
        break
      end
    end

    i = i + 1;
  end

  if found_overlap
    fprintf('The input signals overlap.\n')
  else
    fprintf('The input signals do not overlap.\n')
  end
end


truth.x = zeros(N, numFilters);

switch data_distribution
case DataDistribution.Normal
  for i = 1:numFilters
    truth.x(truth.xSupport(i, :), i) = randn(S, 1);
    truth.x(:, i) = truth.x(:, i) ./ norm(truth.x(:, i));
  end

case DataDistribution.Uniform
  for i = 1:numFilters
    truth.x(truth.xSupport(i, :), i) = rand(S, 1);
    truth.x(:, i) = truth.x(:, i) ./ norm(truth.x(:, i));
  end
end

y = truth.H*truth.x(:);

model.A = kr(model.Psi', model.G.U)';

truth.Zsum = zeros(N, L);
for i = 1:numFilters
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
  truth.Zsum = truth.Zsum + truth.Z{i};
end

end
