% To access the Khatri-Rao product and generate_connected_ER.
addpath('../journal_version/impact_of_rho')

% Number of signals.
P = 1;
% L1 norm minimization.
method = 1;

% Number of nodes in the graph.
N = 20;
% Edge probability in the ER graph.
p = 0.1;
% Build graph, the output is the graph shift operator.
G = generate_connected_ER(N, p);

[V, Lambda] = eig(G);
Vinv = inv(V);
Psi = fliplr(vander(diag(Lambda)));
% Filter order.
L = 3;
Psi = Psi(:,1:L);
% Matrix operator used as a shorthand to compute the output signal y.
M = V*kr(Psi', Vinv')';

% "Sparsity level" of input signals.
S = 2;
% Input signals.
X = randn(S, P);
% Normalize the signals.
X = X./repmat(sqrt(sum(abs(X).^2)), S, 1);
% Pad zeros to generate the full signals.
X = [X; zeros(N-S, P)]; 
x_long = X(:);
x_short = x_long(1:N);

% Generate filter coefficients.
h = randn(L,1); 
% Normalize the filter coefficients.
h = h/norm(h); 

y_short = M*reshape(x_short*h', numel(x_short*h'),1);

% Not used in L1 minimization method.
tau = [];
epsi = [];
dif_epsi = [];

[Z, xEst, hEst, opt] = convex_recovery(y_short, M, P, tau, epsi, dif_epsi, 1);

assert(P == 1)
[Uz, Sz, Vz] = svd(Z');
[hEst sqrt(Sz(1))*Uz(:,1)]
[xEst sqrt(Sz(1))*Vz(:,1)]


