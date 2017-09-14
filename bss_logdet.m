function [Z1_hat, Z2_hat] = bss_logdet(y, A, V, verbose)

if ~exist('verbose', 'var')
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
Theta11_old = eye(N);
Theta12_old = eye(L);
Theta21_old = eye(N);
Theta22_old = eye(L);
% LS solution for initializating w
Z_old = reshape(pinv(B)*y, N, L);
Z1_old = Z_old/2;
Z2_old = Z_old/2;

% Majorization-minimization
while (flag == 1 && iter <= max_iter)
    if verbose
        fprintf('iteration %d\n', iter)
    end

    wx1 = 1./(sqrt(sum(abs(Z1_old).^2, 2)) + epsilon_normx);
    wx2 = 1./(sqrt(sum(abs(Z2_old).^2, 2)) + epsilon_normx);

    cvx_begin quiet
        variable Z1(N, L);
        variable Z2(N, L);
        variable Theta11(N, N) symmetric;
        variable Theta12(L, L) symmetric;
        variable Theta21(N, N) symmetric;
        variable Theta22(L, L) symmetric;

        Z = Z1+Z2;
        pho_Z1 = 5;
        pho_Z2 = 0.5;
        tau_Z1 = 0.5;
        tau_Z2 = 1.3;
        minimize( pho_Z1*(trace((Theta11_old + epsilon_rank*eye(N))\Theta11) + ...
                 trace((Theta12_old + epsilon_rank*eye(L))\Theta12)) + ...
                  pho_Z2*(trace((Theta21_old + epsilon_rank*eye(N))\Theta21) + ...
                 trace((Theta22_old + epsilon_rank*eye(L))\Theta22)) + ...
                 tau_Z1*wx1'*norms(Z1, 2, 2) + tau_Z2*wx2'*norms(Z2, 2, 2) );

        subject to
            [Theta11 Z1; Z1' Theta12] == semidefinite(N+L);
            [Theta21 Z2; Z2' Theta22] == semidefinite(N+L);
            B*Z(:) == y;
    cvx_end

    if verbose
        a = pho_Z1*(trace((Theta11_old + epsilon_rank*eye(N))\Theta11) + ...
                 trace((Theta12_old + epsilon_rank*eye(L))\Theta12));
        b = pho_Z2*(trace((Theta21_old + epsilon_rank*eye(N))\Theta21) + ...
                 trace((Theta22_old + epsilon_rank*eye(L))\Theta22));
        c = tau_Z1*wx1'*norms(Z1, 2, 2);
        d = tau_Z2*wx2'*norms(Z2, 2, 2);

        fprintf('%d %d %d %d\n', a, b, c, d)
    end

    if isempty(strfind(cvx_status, 'Solved'))
        fname = sprintf('failed_problem_bss_logdet_v%s', ...
                        datestr(now, 'ddmmyyyyHHMMSS'));
        warning(sprintf('cvx_status not Solved, saving %s.', fname))
        save(fname)
        Z1_hat = nan(size(Z1));
        Z2_hat = nan(size(Z2));
        return
    end

    difference = norm(Z - Z_old, 'fro')/norm(Z_old, 'fro');

    if ~isempty(strfind(cvx_status, 'Infeasible'))
        % Stop the algorithm.
        if verbose
            fprintf('Infeasible cvx_status.\n')
        end
        Z1_hat = nan(size(Z1));
        Z2_hat = nan(size(Z2));
        return
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
            Z_old = Z;
            Z1_old = Z1;
            Z2_old = Z2;
            Theta11_old = Theta11;
            Theta12_old = Theta12;
            Theta21_old = Theta21;
            Theta22_old = Theta22;
            iter = iter + 1;
        end
    end
end

Z1_hat = Z1;
Z2_hat = Z2;

end
