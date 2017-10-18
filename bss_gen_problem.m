function  [truth, model, y] = bss_gen_problem(params)

if ~exist('params', 'var')
  params = struct;
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
  L = 3;
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
  S = 3;
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
% Because of the way we are coming up with orthogonal vectors,
% which is fixed to three-dimensional vectors.
assert(L == 3)
% Because we are using three-dimensional vectors for the
% filter coefficients and they must be mutually orthogonal.
assert(2 <= numFilters && numFilters <= 3)

switch data_distribution
case DataDistribution.Normal
  truth.h(:, 1) = randn(L, 1);

  truth.h([1 2], 2) = randn(2, 1);

  if numFilters > 2
    truth.h(1, 3) = randn;
  end

case DataDistribution.Uniform
  truth.h(:, 1) = rand(L, 1);

  truth.h([1 2], 2) = rand(2, 1);

  if numFilters > 2
    truth.h(1, 3) = rand;
  end
end

truth.h(3, 2) = - (truth.h([1 2], 1)' * truth.h([1 2], 2)) / truth.h(3, 1);
if numFilters > 2
  truth.h([2 3], 3) = [truth.h([2 3], 1)'; truth.h([2 3], 2)'] \ ...
                      (-truth.h(1, 3) * [truth.h(1,1); truth.h(1,2)]);
end

for i = 1:numFilters
  truth.h(:, i) = truth.h(:, i) / norm(truth.h(:, i));
end

model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1);

% Build filter matrices.
H = zeros(N, N * numFilters);
for i = 1:numFilters
  Hi = truth.h(1, i)*eye(N);
  for l = 1:L-1
    Hi = Hi + truth.h(l+1, i)*model.G.S^l;
  end

  H(:, N*(i-1)+1:N*i) = Hi;
end

% Input.

truth.xSupport = zeros(numFilters, S);
while true
  for i = 1:numFilters
    truth.xSupport(i, :) = randperm(N, S);
  end

  empty_intersection = true;
  for i = 1:numFilters-1
    for j = i+1:numFilters
      if ~isempty(intersect(truth.xSupport(i, :), truth.xSupport(j, :)))
        empty_intersection = false;
      end
    end
  end

  if empty_intersection
    break
  end
end

truth.x = zeros(N, numFilters);

switch data_distribution
case DataDistribution.Normal
  for i = 1:numFilters
    truth.x(truth.xSupport(i, :), i) = randn(S, 1);
  end

case DataDistribution.Uniform
  for i = 1:numFilters
    truth.x(truth.xSupport(i, :), i) = rand(S, 1);
  end
end

% Normalize input signals.
for i = 1:numFilters
  truth.x(:, i) = truth.x(:, i) / norm(truth.x(:, i), 1);
end

for i = 1:numFilters-1
  for j = i+1:numFilters
    assert(truth.x(:, i)' * truth.x(:, j) == 0)
    assert(abs(truth.h(:, i)' * truth.h(:, j)) < 1e-10)
  end
end

y = H*truth.x(:);

model.A = kr(model.Psi', model.G.U)';

truth.Zsum = zeros(N, L);
for i = 1:numFilters
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
  truth.Zsum = truth.Zsum + truth.Z{i};
end

end
