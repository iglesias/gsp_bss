function [Z1_hat, Z2_hat] = sparse_bss_logdet(y, A, V, taux, tauh, varargin)
% SPARSE_BSS_LOGDET: TODO DOC

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
% LS solution for initializating w
Z_old = reshape(pinv(B)*y, N, L);
Theta11_old = eye(N);
Theta12_old = eye(L);
Theta21_old = eye(N);
Theta22_old = eye(L);

% Majorization-minimization
while (flag == 1 && iter <= max_iter)
    fprintf('iteration %d\n', iter)

    wx = 1./(sqrt(sum(abs(Z_old).^2,2)) + epsilon_normx);
    wh = 1./(sqrt(sum(abs(Z_old).^2,1)) + epsilon_normh);

    cvx_begin quiet
        variable Z1(N, L);
        variable Z2(N, L);
        variable Theta11(N, N) symmetric;
        variable Theta12(L, L) symmetric;
        variable Theta21(N, N) symmetric;
        variable Theta22(L, L) symmetric;

        Z = Z1+Z2;
        minimize(trace((Theta11_old + epsilon_rank*eye(N))\Theta11) + ...
                 trace((Theta12_old + epsilon_rank*eye(L))\Theta12) + ...
                 trace((Theta21_old + epsilon_rank*eye(N))\Theta21) + ...
                 trace((Theta22_old + epsilon_rank*eye(L))\Theta22) + ...
                 taux*wx'*norms(Z,2,2) + tauh*norms(Z,2,1)*wh');

        subject to
            [Theta11 Z1; Z1' Theta12] == semidefinite(N+L);
            [Theta21 Z2; Z2' Theta22] == semidefinite(N+L);
            B*Z(:) == y;
    cvx_end

    if isempty(strfind(cvx_status, 'Solved'))
        save(sprintf('failed_problem_sparse_bss_logdet_v%s', ...
                     datestr(now, 'ddmmyyyyHHMMSS')))
    end

    difference = norm(Z - Z_old, 'fro')/norm(Z_old, 'fro');

    if ~isempty(strfind(cvx_status, 'Infeasible'))
        % Stop the algorithm.
        fprintf('Infeasible cvx\n')
        flag = 0;
    else
        if difference < 1e-4
            % Converged.
            fprintf('Convergence reached, cvx_status: %s.\n', cvx_status)
            flag = 0;
        else
            % Did not converge.
            fprintf('Convergence NOT reached, difference=%d.\n', difference)
            Z_old = Z;
            iter = iter + 1;
        end
    end
end

Z1_hat = Z1;
Z2_hat = Z2;

end
