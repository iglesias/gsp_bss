N = 2;
S = [0 1; 1 0];
assert(all(size(S) == [N,N]))

[V, Lambda] = eig(S);
U = eye(N) \ V;

num_filters = 2;
% Number of filters coefficients in each filter.
L = [2 2]';
assert(numel(L) == num_filters)
% Filter coefficients.
h = {};
for i = 1:num_filters
  h{i} = randn(2,1);
end

% Filters.
%TODO generalize.
H1 = h{1}(1) + h{1}(2)*S;
H2 = h{2}(1) + h{2}(2)*S;

Psi1 = repmat(diag(Lambda), 1, L(1)).^[0:L(1)-1];
Psi2 = repmat(diag(Lambda), 1, L(2)).^[0:L(2)-1];

% Input signals.
x1_truth = rand(2,1);
x2_truth = rand(2,1);

% Output signal.
y = H1*x1_truth + H2*x2_truth;
ytilde = [diag(Psi1*h{1}) diag(Psi2*h{2})] * [U*x1_truth; U*x2_truth];

cvx_begin
  variable Z1(2,2)
  variable Z2(2,2)

  Z = Z1+Z2;

  minimize(norm_nuc(Z))
  subject to
    ytilde == kr(Psi1', U')'*vec(Z) 
cvx_end

opt = cvx_optval;
