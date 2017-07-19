function bss(epsilon, taux, tauh, verbose)

if nargin < 2 
  taux = 1e-1;
end

if nargin < 3
  tauh = 1e-1;
end

if nargin < 4
  verbose = false;
end

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
hLP = [1 0]';
hHP = [0 1]';

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

sparse_bss_verbose = false;
[Z1_hat, Z2_hat] = sparse_bss_nuclear(y, A, V, epsilon, taux, tauh, sparse_bss_verbose, knownSupportFlag);
%[Z1_hat, Z2_hat] = sparse_bss_logdet(y, A, V, taux, tauh, sparse_bss_verbose, knownSupportFlag);
%[Z1_hat, Z2_hat] = sparse_bss_logdet_jointsum(y, A, V, 1e-1, 0, sparse_bss_verbose, knownSupportFlag);

Z1 = x1(xSupportToEstimate)*hLP';
Z2 = x2(xSupportToEstimate)*hHP';

save(sprintf('data/epsilon/bss_nuclear_epsilon=%.0d_10_01_%s', epsilon, randomstring(20)))

%[Uz1, Sz1, Vz1] = svd(Z1');
%h1FromZ1 = sqrt(Sz1(1,1))*Uz1(:,1);
%x1FromZ1 = sqrt(Sz1(1,1))*Vz1(:,1);

if verbose
  fprintf('\n\n')
  fprintf('Equality constraint test: %d\n', norm(y - V*A*(Z1_hat(:) + Z2_hat(:))))
  fprintf('norm(Z1)=%d\n', norm(Z1))
  fprintf('norm(Z2)=%d\n', norm(Z2))
  fprintf('norm(Z1_hat)=%d\n', norm(Z1_hat))
  fprintf('norm(Z2_hat)=%d\n', norm(Z2_hat))
  fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - Z1))
  fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - Z2))
  fprintf('norm(Z1_hat-Z2)=%d\n', norm(Z1_hat - Z2))
  fprintf('norm(Z2_hat-Z1)=%d\n', norm(Z2_hat - Z1))
  fprintf('norm([Z1+Z2]-[Z1_hat+Z2_hat])=%d\n', norm([Z1+Z2] - [Z1_hat+Z2_hat]))
  fprintf('norm(Z1-Z2)=%d\n', norm(Z1 - Z2))
  fprintf('norm(Z1_hat-Z2_hat)=%d\n', norm(Z1_hat - Z2_hat))
  fprintf('norm_nuc(Z1)+norm_nuc(Z2)=%d\n', norm_nuc(Z1)+norm_nuc(Z2))
  fprintf('norm_nuc(Z1_hat)=%d\n', norm_nuc(Z1_hat))
  fprintf('norm_nuc(Z2_hat)=%d\n', norm_nuc(Z2_hat))
  fprintf('norm_nuc(Z1_hat)+norm_nuc(Z2_hat)=%d\n', norm_nuc(Z1_hat)+norm_nuc(Z2_hat))
  fprintf('\n\n')

  cminmax = minmax([Z1(:); Z2(:); Z1_hat(:); Z2_hat(:)]);

  subplot(221)
  imagesc(Z1, cminmax), title('Z1')
  subplot(222)
  imagesc(Z2, cminmax), title('Z2')
  subplot(223)
  imagesc(Z1_hat, cminmax), title('Z1\_hat')
  subplot(224)
  imagesc(Z2_hat, cminmax), title('Z2\_hat')
end

end
