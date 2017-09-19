function  [truth, model, y] = bss_gen_problem

% Number of nodes.
N = 50;
% Edge existence probability.
p = 0.1;
% Adjacency matrix.
model.G.W = generate_connected_ER(N, p);
% Graph Laplacian.
model.G.L = diag(sum(model.G.W))-model.G.W;

assert(issymmetric(model.G.W))
assert(issymmetric(model.G.L))

%figure
%matlabGraph = graph(G.W);
%plot(matlabGraph)

[model.G.V, model.G.D] = eig(model.G.L);
model.G.U = inv(model.G.V);
model.G.lambda = diag(model.G.D);

numFilterCoeffs = 2;
data_distibution = DataDistribution.Uniform;

switch data_distibution
case DataDistribution.Normal
  truth.h1 = randn(2, 1);
  truth.h2 = randn(2, 1);

  truth.h1 = truth.h1/norm(truth.h1);
  truth.h2 = truth.h2/norm(truth.h2);
case DataDistribution.Uniform
  truth.h1 = rand(2, 1);
  truth.h2 = rand(2, 1);
end

Psi = repmat(model.G.lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);

% Build filter matrices.
H1 = truth.h1(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H1 = H1 + truth.h1(l+1)*model.G.L^l;
end

H2 = truth.h2(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H2 = H2 + truth.h2(l+1)*model.G.L^l;
end

H = [H1 H2];

% Input.
% Number of non-zero input nodes.
S = 5;

truth.x1Support = randperm(N, S);
truth.x2Support = randperm(N, S);

if isempty(intersect(truth.x1Support, truth.x2Support))
  fprintf('Intersection of the inputs'' supports empty.\n')
else
  fprintf('Intersection of the inputs'' supports non-empty.\n')
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

model.A = kr(Psi', model.G.U')';
xSupportToEstimate = 1:N;

truth.Z1 = truth.x1(xSupportToEstimate)*truth.h1';
truth.Z2 = truth.x2(xSupportToEstimate)*truth.h2';

end
