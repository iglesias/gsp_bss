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
lambda = diag(D);

showFlag = 0;
numFilterCoeffs = 10;
[hLP, hHP] = gsp_design_filters(lambda, numFilterCoeffs, showFlag);

% Build filter matrix.
H = hLP(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H = H + hLP(l+1)*G.L^l;
end

% Input.
x = zeros(N, 1);
% Number of non-zero input nodes.
S = 20;
xSupport = randperm(N, S);
x(xSupport) = rand(S, 1);

y = H*x;
xest = H\y;
