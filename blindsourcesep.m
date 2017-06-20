clearvars

% Number of nodes.
N = 20;
% Edge existence probability.
p = 0.1;
% Adjacency matrix.
G.W = generate_connected_ER(N, p);
% Graph Laplacian.
G.L = diag(sum(G.W))-G.W;

assert(issymmetric(G.W))
assert(issymmetric(G.L))

%figure
%matlabGraph = graph(G.W);
%plot(matlabGraph)

[V, D] = eig(G.L);
U = inv(V);
lambda = diag(D);

showFlag = 0;
numFilterCoeffs = 3;
[hLP, hHP] = gsp_design_filters(lambda, numFilterCoeffs, showFlag);

Psi = repmat(lambda, 1, numFilterCoeffs).^[0:numFilterCoeffs-1];

% Build filter matrices.
H1 = hLP(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H1 = H1 + hLP(l+1)*G.L^l;
end

H2 = hHP(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H2 = H2 + hHP(l+1)*G.L^l;
end

H = [H1 H2];

I = eye(N);

% Input.
% Number of non-zero input nodes.
S = 2;

xSupport = randperm(N, S);

x1 = zeros(N, 1);
x1(xSupport) = rand(S, 1);

x2 = zeros(N, 1);
x2(xSupport) = rand(S, 1);

x = [x1; x2];
y = H*x;

Ubar = U(:, xSupport);
A = kr(Psi', Ubar')';

Z1 = x1(xSupport)*hLP';
Z2 = x2(xSupport)*hHP';
knownSupportFlag = 1;
[Z1est, Z2est] = sparse_bss_nuclear(y, A, V, 1e-1, 0, knownSupportFlag, S);

[Z1 Z2]
[Z1est Z2est]

[Z1+Z2]
[Z1est+Z2est]
