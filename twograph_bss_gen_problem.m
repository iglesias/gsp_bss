numGraphs = 2;
% Number of nodes.
N = 50;
% Edge existence probability.
p = 0.1;
% Adjacency matrices.
for i = 1:numGraphs
    G(i).W = generate_connected_ER(N, p);
    G(i).L = diag(sum(G(i).W))-G(i).W;
    [G(i).V, G(i).D] = eig(G(i).L);
    G(i).U = inv(G(i).V);
    G(i).lambda = diag(G(i).D);
end

numFilterCoeffs = 2;
h1 = rand(2, 1);
h2 = rand(2, 1);

for i = 1:numGraphs
  Psi{i} = repmat(G(i).lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);
end

% Build filter matrices.
H1 = h1(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H1 = H1 + h1(l+1)*G(1).L^l;
end

H2 = h2(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H2 = H2 + h2(l+1)*G(2).L^l;
end

H = [H1 H2];

% Input.
% Number of non-zero input nodes.
S = 5;

global x1Support x2Support
x1Support = randperm(N, S);
x2Support = randperm(N, S);

fprintf('Intersection of the supports of the inputs:\n')
intersect(x1Support, x2Support)

x1 = zeros(N, 1);
x1(x1Support) = rand(S, 1);

x2 = zeros(N, 1);
x2(x2Support) = rand(S, 1);

x = [x1; x2];
y = H*x;

A1 = kr(Psi{1}', G(1).U')';
A2 = kr(Psi{2}', G(2).U')';

Z1 = x1*h1';
Z2 = x2*h2';
