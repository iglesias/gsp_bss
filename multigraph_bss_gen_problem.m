function [truth, model, y] = multigraph_bss_gen_problem(params)

data_distribution = DataDistribution.Normal;
shift_operator = ShiftOperator.Adjacency;

if ~exist('params', 'var')
  params = struct;
end

if isfield(params, 'numGraphs')
  R = params.numGraphs;
else
  R = 2;
end

% Number of filter coefficients.
if isfield(params, 'L')
  L = params.L;
else
  L = 3;
end

if isfield(params, 'filterDistribution')
  filter_distribution = params.filterDistribution;
else
  filter_distribution = DataDistribution.Uniform;
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

% Edge existence probability.
p = 0.1;
for i = 1:R
  % Adjacency matrix.
  model.G(i).W = generate_connected_ER(N, p);
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

  [model.G(i).V, Lambda] = eig(model.G(i).S);
  model.G(i).U = inv(model.G(i).V);
  model.G(i).lambda = diag(Lambda);
  % disp(minmax(model.G(i).lambda'))
end

model.V = reshape([model.G.V], [N, N, R]);

% Filter coefficients.
truth.h = zeros(L, R);

switch filter_distribution
case DataDistribution.Normal
  truth.h = randn(L, R);
case DataDistribution.Uniform
  truth.h = rand(L, R);
case DataDistribution.HeatKernel
  truth.h = exp(-(0.3*rand(R, 1)+0.7)*(L:-1:1))';
end

% Normalize filter taps.
truth.h = truth.h ./ repmat(norms(truth.h, 2, 1), L, 1);

for i = 1:R
  model.Psi{i} = repmat(model.G(i).lambda, 1, L).^repmat([0:L-1], N, 1);
  % disp(minmax(vec(model.Psi{i})'))
end

% Build filter matrices.
H = zeros(N, N * R);
for i = 1:R
  Hi = truth.h(1, i)*eye(N);
  for l = 1:L-1
    Hi = Hi + truth.h(l+1, i)*model.G(i).S^l;
  end

  H(:, N*(i-1)+1:N*i) = Hi;
end

% Input signals.

truth.xSupport = zeros(R, S);
for i = 1:R
  truth.xSupport(i, :) = randperm(N, S);
end

truth.x = zeros(N, R);

switch data_distribution
case DataDistribution.Normal
  for i = 1:R
    truth.x(truth.xSupport(i, :), i) = randn(S, 1);
  end

case DataDistribution.Uniform
  % TODO vectorize.
  for i = 1:R
    truth.x(truth.xSupport(i, :), i) = rand(S, 1);
  end
end

% Normalize input signals.
% TODO vectorize.
for i = 1:R
  truth.x(:, i) = truth.x(:, i) ./ norm(truth.x(:, i));
end

y = H*truth.x(:);

for i = 1:R
  model.A{i} = kr(model.Psi{i}', model.G(i).U')';
  % disp(minmax(vec(model.A{i})'))
end

% TODO vectorize.
for i = 1:R
  truth.Z{i} = truth.x(:, i)*truth.h(:, i)';
end

end
