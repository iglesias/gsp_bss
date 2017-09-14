function [Z1_hat, Z2_hat] = bss_logdet_jointsum(y, A, V, taux, tauh, verbose, varargin)
% BSS_LOGDET_JOINTSUM: TODO DOC

if nargin < 6
    verbose = false;
end

%% Parameter definition.

epsilon_rank = 5e-2;
epsilon_normx = 1e-3;
epsilon_normh = 1e-3;
max_iter = 50;

%% Blind source separation using rank mininimization: log-det surrogate

assert(size(A, 1) == length(y))
N = size(A, 1);
% It is assumed that the filters have the same order.
L = size(A, 2)/N;
B = V*A;

% Initialization
flag = 1;
iter = 1;
Theta1_old = eye(N);
Theta2_old = eye(L);
% LS solution for initializating w
Z_old = reshape(pinv(B)*y, N, L);

% Majorization-minimization
while (flag == 1 && iter <= max_iter)
    if verbose
        fprintf('iteration %d\n', iter)
    end

    wx = 1./(sqrt(sum(abs(Z_old).^2, 2)) + epsilon_normx);
    wh = 1./(sqrt(sum(abs(Z_old).^2, 1)) + epsilon_normh);

    cvx_begin quiet
        variable Z1(N, L);
        variable Z2(N, L);
        variable Theta1(N, N) symmetric;
        variable Theta2(L, L) symmetric;

        Z = Z1+Z2;
        minimize(trace((Theta1_old + epsilon_rank*eye(N))\Theta1) + ...
                 trace((Theta2_old + epsilon_rank*eye(L))\Theta2) + ...
                 taux*wx'*norms(Z,2,2) + tauh*norms(Z,2,1)*wh');

        subject to
            [Theta1 Z; Z' Theta2] == semidefinite(N+L);
            B*Z(:) == y;
    cvx_end

    if isempty(strfind(cvx_status, 'Solved'))
        fname = sprintf('failed_problem_bss_logdet_jointsum_v%s', ...
                        datestr(now, 'ddmmyyyyHHMMSS'));
        warning(sprintf('cvx_status not Solved, saving %s.', fname))
        save(fname)
    end

    difference = norm(Z - Z_old, 'fro')/norm(Z_old, 'fro');

    if ~isempty(strfind(cvx_status, 'Infeasible'))
        % Stop the algorithm.
        if verbose
            fprintf('Infeasible cvx_status.\n')
        end
        flag = 0;
    else
        if difference < 1e-4
            % Converged.
            if verbose
                fprintf('Convergence reached, cvx_status: %s.\n', cvx_status)
            end
            flag = 0;
        else
            % Did not converge.
            if verbose
                fprintf('Convergence NOT reached, difference=%d.\n', difference)
            end
            Theta1_old = Theta1;
            Theta2_old = Theta2;
            Z_old = Z;
            iter = iter + 1;
        end
    end
end

Z1_hat = Z1;
Z2_hat = Z2;

end
