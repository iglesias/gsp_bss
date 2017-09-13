function [Z1_hat, Z2_hat] = twograph_bss_nuclear_direct(y, A1, V1, A2, V2, verbose)

if ~exist('verbose', 'var')
  verbose = false;
end

assert(all(size(A1) == size(A2)))
assert(size(A1, 1) == length(y))
N = size(A1, 1);
% It is assumed that the filters have the same order.
L = size(A1, 2)/N;

cvx_begin quiet
    variable Z1(N, L);
    variable Z2(N, L);

    pho = 1;
    tau = 0.5;
    minimize(pho*(norm_nuc(Z1) + norm_nuc(Z2)) + tau*(sum(norms(Z1, 1, 2)) + sum(norms(Z2, 1, 2))));

    subject to
        y == V1*A1*vec(Z1) + V2*A2*vec(Z2);
cvx_end

if verbose
    fprintf('%d %d\n', pho*(norm_nuc(Z1) + norm_nuc(Z2)), tau*(sum(norms(Z1, 1, 2)) + sum(norms(Z2, 1, 2))))
end

Z1_hat = Z1;
Z2_hat = Z2;

end
