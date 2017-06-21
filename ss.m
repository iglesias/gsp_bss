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
S = 5;

x1 = zeros(N, 1);
x1Support = randperm(N, S);
x1(x1Support) = rand(S, 1);

x2 = zeros(N, 1);
x2Support = randperm(N, S);
x2(x2Support) = rand(S, 1);

x = [x1; x2];
y = H*x;

xest = knownh_bci(y, H);
norm(xest - x)
