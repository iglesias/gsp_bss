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

shift = G.L;

[V, D] = eig(shift);
U = inv(V);
lambda = diag(D);

showFlag = 0;
numFilterCoeffs = 10;
[hLP, hHP] = gsp_design_filters(lambda, numFilterCoeffs, showFlag);

% Build filter matrices.
H1 = hLP(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H1 = H1 + hLP(l+1)*shift^l;
end

H2 = hHP(1)*eye(N);
for l = 1:numFilterCoeffs-1
  H2 = H2 + hHP(l+1)*shift^l;
end

S = 2;
x = zeros(N, 1);
xSupport = sort(randperm(N, S));
x(xSupport) = rand(S, 1);
x = x / norm(x);

xtilde = U*x;
y = H2*x;
ytilde = U*y;

figure(1)
subplot(211)
stem(x)
subplot(212)
stem(xtilde)
axis([1 N min([xtilde; ytilde]) max([xtilde; ytilde])])

figure(2)
subplot(211)
stem(y)
subplot(212)
stem(ytilde)
axis([1 N min([xtilde; ytilde]) max([xtilde; ytilde])])
