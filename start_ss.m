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

figure
matlabGraph = graph(G.W);
plot(matlabGraph)

[V, D] = eig(G.L);
lambda = diag(D);

showFlag = 1;
gsp_design_filters(lambda, 10, showFlag);
