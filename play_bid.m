N = 2;
S = [0 1; 1 0];
assert(all(size(S) == [N,N]))

[V, Lambda] = eig(S);
% Signal gft.
U = eye(N) \ V;

% Number of filter coefficients.
L = 2;
h = randn(2,1);

% Filter.
H = h(1) + h(2)*S;

% Filter gft.
Psi = repmat(diag(Lambda), 1, L).^[0:L-1];

% Input signal.
x_truth = rand(2,1);

% Output signal.
y = H*x_truth;
% In frequency.
ytilde = diag(Psi*h) * (U*x_truth);

% Build matrix that will be useful for the contraint with ytilde.
A = kr(Psi', U')';

cvx_begin
  variable Z(2,2)

  minimize(norm_nuc(Z))
  subject to
    ytilde == A*vec(Z)
cvx_end

opt = cvx_optval;
