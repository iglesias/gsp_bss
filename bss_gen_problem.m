% Number of nodes.
N = 50;
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

numFilterCoeffs = 2;
h1 = rand(2, 1);
h2 = rand(2, 1);

Psi = repmat(lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);

% Build filter matrices.
H1 = h1(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H1 = H1 + h1(l+1)*G.L^l;
end

H2 = h2(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H2 = H2 + h2(l+1)*G.L^l;
end

H = [H1 H2];

I = eye(N);

% Input.
% Number of non-zero input nodes.
S = 5;

x1Support = randperm(N, S);
x2Support = randperm(N, S);

%intersect(x1Support, x2Support)

x1 = zeros(N, 1);
x1(x1Support) = rand(S, 1);

x2 = zeros(N, 1);
x2(x2Support) = rand(S, 1);

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

Z1 = x1(xSupportToEstimate)*h1';
Z2 = x2(xSupportToEstimate)*h2';
