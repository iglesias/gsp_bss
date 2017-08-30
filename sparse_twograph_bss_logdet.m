function [Z1_hat, Z2_hat] = sparse_twograph_bss_logdet(y, A1, V1, A2, V2, verbose)

if ~exist('verbose', 'var')
    verbose = false;
end

%% Parameter definition.

epsilon_rank = 5e-2;
epsilon_normx = 1e-3;
max_iter = 50;

%% Blind source separation using rank mininimization: log-det surrogate
assert(all(size(A1) == size(A2)))
assert(size(A1, 1) == length(y))
N = size(A1, 1);
% It is assumed that the filters have the same order.
L = size(A1, 2)/N;

% Initialization
flag = 1;
iter = 1;
Theta11_old = eye(N);
Theta12_old = eye(L);
Theta21_old = eye(N);
Theta22_old = eye(L);
%FIXME does this initialization make sense?
% Initial solution for initializing w
global x1Support x2Support
%x1 = zeros(N, 1);
%x1(x1Support) = rand(length(x1Support), 1);
%x1(randperm(N, 5)) = rand(5, 1);
x1 = rand(N, 1);
%x2 = zeros(N, 1);
%x2(x2Support) = rand(length(x2Support), 1);
%x2(randperm(N, 5)) = rand(5, 1);
x2 = rand(N, 1);
Z1_old = x1*rand(1, L);
Z2_old = x2*rand(1, L);

% Majorization-minimization
pho_Z1 = 1;
pho_Z2 = 1;
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

        pho_Z1 = 10;
        pho_Z2 = 10;
        minimize(pho_Z1*(trace((Theta11_old + epsilon_rank*eye(N))\Theta11) + ...
                 trace((Theta12_old + epsilon_rank*eye(L))\Theta12)) + ...
                 pho_Z2*(trace((Theta21_old + epsilon_rank*eye(N))\Theta21) + ...
                 trace((Theta22_old + epsilon_rank*eye(L))\Theta22)) + ...
                 wx1'*norms(Z1, 2, 2) + wx2'*norms(Z2, 2, 2));

        subject to
            [Theta11 Z1; Z1' Theta12] == semidefinite(N+L);
            [Theta21 Z2; Z2' Theta22] == semidefinite(N+L);
            y == V1*A1*vec(Z1) + V2*A2*vec(Z2);
    cvx_end

    a = pho_Z1*(trace((Theta11_old + epsilon_rank*eye(N))\Theta11) + ...
                trace((Theta12_old + epsilon_rank*eye(L))\Theta12)) + ...
        pho_Z2*(trace((Theta21_old + epsilon_rank*eye(N))\Theta21) + ...
                trace((Theta22_old + epsilon_rank*eye(L))\Theta22));
    b = wx1'*norms(Z1, 2, 2) + wx2'*norms(Z2, 2, 2);
    fprintf('%d %d\n', a, b)

    if isempty(strfind(cvx_status, 'Solved'))
        fname = sprintf('failed_problem_sparse_twograph_bss_logdet_v%s', ...
                        datestr(now, 'ddmmyyyyHHMMSS'));
        warning(sprintf('cvx_status not Solved, saving %s.', fname))
        save(fname)
        Z1_hat = nan(size(Z1));
        Z2_hat = nan(size(Z2));
        return
    end

    difference = norm(Z1 - Z1_old, 'fro')/norm(Z1_old, 'fro') + ...
                 norm(Z2 - Z2_old, 'fro')/norm(Z2_old, 'fro');

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
