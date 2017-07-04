clearvars

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

%showFlag = 0;
%numFilterCoeffs = 15;
%[hLP, hHP] = gsp_design_filters(lambda, numFilterCoeffs, showFlag);

numFilterCoeffs = 2;
hLP = [0 1]';
hHP = [-1 0]';

Psi = repmat(lambda, 1, numFilterCoeffs).^repmat([0:numFilterCoeffs-1], N, 1);

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

x1Support = randperm(N, S);
x2Support = randperm(N, S);

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

[Z1est, Z2est] = sparse_bss_nuclear(y, A, V, 1e-1, 0, knownSupportFlag);

Z1 = x1(xSupportToEstimate)*hLP';
Z2 = x2(xSupportToEstimate)*hHP';

[Uz1, Sz1, Vz1] = svd(Z1');
h1FromZ1 = sqrt(Sz1(1,1))*Uz1(:,1);
x1FromZ1 = sqrt(Sz1(1,1))*Vz1(:,1);

fprintf('\n\n')
fprintf('norm(Z1est-Z1)=%d\n', norm(Z1est - Z1))
fprintf('norm(Z2est-Z2)=%d\n', norm(Z2est - Z2))
fprintf('norm([Z1+Z2]-[Z1est+Z2est])=%d\n', norm([Z1+Z2] - [Z1est+Z2est]))
fprintf('norm(Z1-Z2)=%d\n', norm(Z1 - Z2))
fprintf('norm(Z1est-Z2est)=%d\n', norm(Z1est - Z2est))
fprintf('norm_nuc(Z1)+norm_nuc(Z2)=%d\n', norm_nuc(Z1)+norm_nuc(Z2))
fprintf('norm_nuc(Z1est)+norm_nuc(Z2est)=%d\n', norm_nuc(Z1est)+norm_nuc(Z2est))
fprintf('\n\n')

subplot(221)
imagesc(Z1), title('Z1')
subplot(222)
imagesc(Z2), title('Z2')
subplot(223)
imagesc(Z1est), title('Z1est')
subplot(224)
imagesc(Z2est), title('Z2est')
