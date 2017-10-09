function  [truth, model, y] = brain_bss_gen_problem(brain_idxs)

assert(length(brain_idxs) > 0)

numGraphs = length(brain_idxs);
numFilters = numGraphs;
% Number of non-zero input nodes.
S = 1;

data_distribution = DataDistribution.Uniform;
shift_operator = ShiftOperator.Adjacency;

load('data/brain_data_66', 'CC')
N = 66;

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
% Order of the filters (number of filter coefficients).
L = 3;
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

for i = 1:numGraphs
  model.Psi{i} = repmat(model.G(i).lambda, 1, L).^repmat([0:L-1], N, 1);
end

% Build filter matrices.
H = zeros(N, N * numFilters);
for i = 1:numFilters
  Hi = truth.h(1, i)*eye(N);
  for l = 1:L-1
    Hi = Hi + truth.h(l+1, i)*model.G(i).S^l;
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

for i = 1:numGraphs
  model.A{i} = kr(model.Psi{i}', model.G(i).U)';
end

for i = 1:numGraphs
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
end

end
