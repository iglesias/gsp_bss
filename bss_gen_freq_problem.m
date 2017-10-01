function  [truth, model, y_tilde] = bss_gen_freq_problem

numFilters = 2;
numFilterCoeffs = 2;
data_distibution = DataDistribution.Uniform;
% Number of non-zero input nodes.
S = 1;

% Number of nodes.
N = 10;

% Edge existence probability.
p = 0.1;
% Adjacency matrix.
model.G.W = generate_connected_ER(N, p);
% Graph Laplacian.
model.G.L = diag(sum(model.G.W))-model.G.W;

assert(issymmetric(model.G.W))
assert(issymmetric(model.G.L))

[model.G.V, model.G.D] = eig(model.G.L);
model.G.U = inv(model.G.V);
model.G.lambda = diag(model.G.D);

model.Psi = repmat(model.G.lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);

assert(mod(N, numFilters) == 0)
h1_tilde = [rand(2, 1); zeros(N-2, 1)];
h2_tilde = [zeros(N-2, 1); rand(2, 1)];

% Input.

while true
  truth.x1Support = randperm(N, S);
  truth.x2Support = randperm(N, S);
  if isempty(intersect(truth.x1Support, truth.x2Support))
    break
  end
end

truth.x1 = zeros(N, 1);
truth.x2 = zeros(N, 1);

switch data_distibution
case DataDistribution.Normal
  truth.x1(truth.x1Support) = randn(S, 1);
  truth.x2(truth.x2Support) = randn(S, 1);

case DataDistribution.Uniform
  truth.x1(truth.x1Support) = rand(S, 1);
  truth.x2(truth.x2Support) = rand(S, 1);
end

truth.x1 = truth.x1 / norm(truth.x1, 1);
truth.x2 = truth.x2 / norm(truth.x2, 1);

x1_tilde = model.G.U * truth.x1;
x2_tilde = model.G.U * truth.x2;

y_tilde = h1_tilde .* x1_tilde + h2_tilde .* x2_tilde;

model.A = kr(eye(N), model.G.U')';

truth.ZPsi1 = truth.x1*h1_tilde';
truth.ZPsi2 = truth.x2*h2_tilde';

[UZPsi, SZPsi, VZPsi] = svd(truth.ZPsi1 + truth.ZPsi2);
norm(SZPsi(1,1)*UZPsi(:,1)*VZPsi(:,1)' - (truth.x1 * h1_tilde'))
norm(SZPsi(2,2)*UZPsi(:,2)*VZPsi(:,2)' - (truth.x2 * h2_tilde'))
norm(SZPsi(1,1)*UZPsi(:,1)*VZPsi(:,1)' - (truth.x2 * h2_tilde'))
norm(SZPsi(2,2)*UZPsi(:,2)*VZPsi(:,2)' - (truth.x1 * h1_tilde'))

end
