clearvars

% Number of nodes.
N = 100;
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
numFilterCoeffs = 15;
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
S = 10;

xSupport = randperm(N, S);

x1 = zeros(N, 1);
x1(xSupport) = rand(S, 1);

x2 = zeros(N, 1);
x2(xSupport) = rand(S, 1);

x = [x1; x2];
y = H*x;

knownSupportFlag = false;

if knownSupportFlag
  % Known support equal for all inputs
  % (see (3) in CAMSAP_BlindId for the expression below in the single-input case).
  Ubar = U(:, xSupport);
  A = kr(Psi', Ubar')';
  xSupportToEstimate = xSupport;
else
  A = kr(Psi', U')';
  xSupportToEstimate = 1:N;
end

[Z1est, Z2est] = sparse_bss_nuclear(y, A, V, 1e-1, 0, knownSupportFlag);

Z1 = x1(xSupportToEstimate)*hLP';
Z2 = x2(xSupportToEstimate)*hHP';
norm(Z1est - Z1)
norm(Z2est - Z2)
norm([Z1+Z2] - [Z1est+Z2est])
norm(Z1est - Z2est)
