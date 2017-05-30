%%% Graph.

% Number of nodes.
N = 20;
% Edge probability in the Erdos-Renyi graph.
p = 0.1;
% Build graph, the output is the graph-shift operator.
S = generate_connected_ER(N, p);
assert(all(size(S) == [N,N]))

[V, Lambda] = eig(S);
U = inv(V);

num_filters = 2;
% Number of filters coefficients in each filter.
L = [2 2]';
assert(numel(L) == num_filters)
% Filter coefficients.
h = {};
for i = 1:num_filters
  h{i} = randn(2,1);
end

%%% Filters.

%TODO generalize to Li filter coefficients.
assert(all(L == 2))
%TODO generalize to arbitraty number of filters.
assert(num_filters == 2)
H1 = h{1}(1) + h{1}(2)*S;
H2 = h{2}(1) + h{2}(2)*S;
Psi1 = repmat(diag(Lambda), 1, L(1)).^[0:L(1)-1];
Psi2 = repmat(diag(Lambda), 1, L(2)).^[0:L(2)-1];

%%% Input.

x1_truth = zeros(N, 1);
x2_truth = zeros(N, 1);

% nnz elements.
S = 3;

x1_support = randperm(N, S);
x2_support = randperm(N, S);

x1_truth(x1_support) = randn(S, 1);
x2_truth(x2_support) = randn(S, 1);

x1_truth = x1_truth/norm(x1_truth);
x2_truth = x2_truth/norm(x2_truth);

%%% Output.

y = H1*x1_truth + H2*x2_truth;
% In frequency.
ytilde = [diag(Psi1*h{1}) diag(Psi2*h{2})] * [U*x1_truth; U*x2_truth];

%%% Reconstruction.

% Build matrix that will be useful for the constraint with ytilde.
% If the degree of the filters is not equal, then their Psi are different
% and the expression in the whiteboard should change.
assert(L(1) == L(2))
A = kr(Psi1', U')';

cvx_begin
  variable Z1(N, L(1))
  variable Z2(N, L(2))

  Z = Z1+Z2;

  minimize(norm_nuc(Z))
  subject to
    ytilde == A*vec(Z)
cvx_end

opt = cvx_optval;
