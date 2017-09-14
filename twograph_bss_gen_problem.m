function [truth, model, y] = twograph_bss_gen_problem

numGraphs = 2;
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

numFilterCoeffs = 3;
data_distibution = DataDistribution.Uniform;

switch data_distibution
case DataDistribution.Normal
  truth.h1 = randn(numFilterCoeffs, 1);
  truth.h2 = randn(numFilterCoeffs, 1);

  truth.h1 = truth.h1/norm(truth.h1);
  truth.h2 = truth.h2/norm(truth.h2);
case DataDistribution.Uniform
  truth.h1 = rand(numFilterCoeffs, 1);
  truth.h2 = rand(numFilterCoeffs, 1);
end

for i = 1:numGraphs
  Psi{i} = repmat(model.G(i).lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);
end

% Build filter matrices.
H1 = truth.h1(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H1 = H1 + truth.h1(l+1)*model.G(1).L^l;
end

H2 = truth.h2(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H2 = H2 + truth.h2(l+1)*model.G(2).L^l;
end

H = [H1 H2];

% Input.
% Number of non-zero input nodes.
S = 1;

truth.x1Support = randperm(N, S);
truth.x2Support = randperm(N, S);

if isempty(intersect(truth.x1Support, truth.x2Support))
  fprintf('Intersection of the supports of the inputs empty.\n')
else
  fprintf('Intersection of the supports of the inputs non-empty.\n')
end

truth.x1 = zeros(N, 1);
truth.x2 = zeros(N, 1);

switch data_distibution
case DataDistribution.Normal
  truth.x1(truth.x1Support) = randn(S, 1);
  truth.x2(truth.x2Support) = randn(S, 1);

  truth.x1 = truth.x1/norm(truth.x1);
  truth.x2 = truth.x2/norm(truth.x2);

case DataDistribution.Uniform
  truth.x1(truth.x1Support) = rand(S, 1);
  truth.x2(truth.x2Support) = rand(S, 1);
end


x = [truth.x1; truth.x2];
y = H*x;

model.A1 = kr(Psi{1}', model.G(1).U')';
model.A2 = kr(Psi{2}', model.G(2).U')';

truth.Z1 = truth.x1*truth.h1';
truth.Z2 = truth.x2*truth.h2';

end
