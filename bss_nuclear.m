function [Z1_hat, Z2_hat] = bss_nuclear(y, A, V, verbose)
% BSS_NUCLEAR: Obtains the estimate of the filter and the elements
% of a graph signal using the nuclear norm surrogate, assuming a sparse model for the signal
%
%           [Z1_hat, Z2_hat] = bss_nuclear(y, A, V, verbose)
%           y = filtered graph signal (might be a subset)
%           A = Khatri-Rao product matrix
%           V = rows of the inverse Graph Fourier transform which correspond to the observations
%           verbose = set to true/false to switch on/off console output (default false)
%
%           Zi_hat = estimate of the rank-1 matrices xi*hi'

%% Parameter definition

epsilon_normx = 1e-3;
epsilon_normh = 1e-3;
maxiter = 50;

%% Check input parameters

if ~exist('verbose', 'var')
    verbose = false;
end

%% Blind source separation using rank minimization: nuclear norm surrogate

assert(size(A, 1) == length(y)) % this might change for sampled observations
N = size(A, 1);
% It is assumed that the filters have the same order.
L = size(A, 2)/N;
B = V*A;

% Initialization
flag = 1;
iter = 1;
% LS solution for initializing w
Z_old = reshape(pinv(B)*y, N, L);
Z1_old = Z_old/2;
Z2_old = Z_old/2;

% Majorization-minimization
while (flag == 1 && iter <= maxiter)
    if verbose
        fprintf('iteration %d\n', iter)
    end

    wx1 = 1./(sqrt(sum(abs(Z1_old).^2, 2)) + epsilon_normx);
    wx2 = 1./(sqrt(sum(abs(Z2_old).^2, 2)) + epsilon_normx);

    cvx_begin quiet
        variable Z1(N, L);
        variable Z2(N, L);

        Z = Z1+Z2;
        pho_Z1 = 5;
        tau_Z2 = 1.3;
        minimize( pho_Z1*norm_nuc(Z1) + norm_nuc(Z2) + ...
                  wx1'*norms(Z1, 2, 2) + tau_Z2*wx2'*norms(Z2, 2, 2) );

        subject to
            B*Z(:) == y;
    cvx_end

%    fprintf('%d %d %d %d\n', pho_Z1*norm_nuc(Z1), norm_nuc(Z2), ...
%                             wx1'*norms(Z1, 2, 2), ...
%                             tau_Z2*wx2'*norms(Z2, 2, 2))

    if isempty(strfind(cvx_status, 'Solved'))
        fname = sprintf('failed_problem_bss_nuclear_v%s', ...
                        datestr(now, 'ddmmyyyyHHMMSS'));
        warning(sprintf('cvx_status not Solved, saving %s.', fname))
        save(fname)
        Z1_hat = nan(size(Z1));
        Z2_hat = nan(size(Z2));
        return
    end

    difference = norm(Z - Z_old,'fro')/norm(Z_old,'fro');

    if ~isempty(strfind(cvx_status, 'Infeasible'))
        % Stop the algorithm
        if verbose
            fprintf('Infeasible cvx\n')
        end
        Z1_hat = nan(size(Z1));
        Z2_hat = nan(size(Z2));
        return
    else
        if difference < 1e-4
            % Converged
            if verbose
                fprintf('Convergence reached, difference=%d, cvx_status: %s.\n', ...
                        difference, cvx_status)
            end
            flag = 0;
        else
            % Did not converge
            if verbose
                fprintf('Covergence NOT reached, difference=%d.\n', difference)
            end
            Z_old = Z;
            Z1_old = Z1;
            Z2_old = Z2;
            iter = iter + 1;
        end
    end

end

Z1_hat = Z1;
Z2_hat = Z2;

end
