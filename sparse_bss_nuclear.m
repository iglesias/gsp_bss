function [Z1_hat, Z2_hat] = sparse_bss_nuclear(y, A, V, taux, tauh, verbose, varargin)
% SPARSE_BSS_NUCLEAR: Obtains the estimate of the filter and the elements 
% of a graph signal using the nuclear norm surrogate, assuming a sparse model for the signal
%
%
%           [Z1_hat, Z2_hat] = sparse_bss_nuclear(y, M, L, taux, tauh, verbose, x_known, known_indexes)
%           y = filtered graph signal (might be a subset)
%           A = Khatri-Rao product matrix
%           V = rows of the inverse Graph Fourier transform which correspond to the observations
%           taux = regularization parameter for x
%           tauh = regularization parameter for h
%           known_support = whether the support of the input signals (assumed the same) is known
%           S = cardinality of the input signals' support
%           x_known = known values of the graph signal x. Could be an empty matrix or optional
%           known_indexes = indexes of the known values. Could be an empty matrix or optional
%       
%           Zi_hat = estimate of the rank-1 matrices xi*hi'

warning ('off','MATLAB:nargchk:deprecated')

%% Parameter definition

epsilon_normx = 1e-3;
epsilon_normh = 1e-3;
maxiter = 50;

%% Check input parameters

if nargin < 6
    verbose = false;
end

if nargin < 7
    known_support = 0;
    S = [];
else 
    assert(nargin==7 || nargin==9)
    known_support = varargin{1};
    if known_support
        S = varargin{2};
    end
end

if nargin < 9
    known_indexes = [];
    x_known = [];
else
    assert(nargin==9)
    x_known = varargin{3};
    known_indexes = varargin{4};
end

%% Blind source separation using rank minimization: nuclear norm surrogate

assert(size(A, 1) == length(y)) % this might change for sampled observations
N = size(A, 1);
% It is assumed that the filters have the same order.
if known_support
    L = size(A, 2)/S;
else
    L = size(A, 2)/N;
end
B = V*A;

% Length of the known signal values
N_known = length(known_indexes);

% Initialization
flag = 1;
iter = 1;
% LS solution for initializing w
if known_support
    Z_old = reshape(pinv(B)*y,S,L);
else
    Z_old = reshape(pinv(B)*y,N,L);
    Z1_old = Z_old/2;
    Z2_old = Z_old/2;
end

% Majorization-minimization
while (flag == 1 && iter <= maxiter)
    if verbose
        fprintf('iteration %d\n', iter)
    end

    wx = 1./(sqrt(sum(abs(Z_old).^2,2)) + epsilon_normx);
    wx(known_indexes) = 0;
    %wh = 1./(sqrt(sum(abs(Z_old).^2,1)) + epsilon_normh);
    wh1 = 1./(sqrt(sum(abs(Z1_old).^2,1)) + epsilon_normh);
    wh2 = 1./(sqrt(sum(abs(Z2_old).^2,1)) + epsilon_normh);
    
    cvx_begin quiet
        if known_support
            variable Z1(S,L);
            variable Z2(S,L);
        else
            variable Z1(N,L);
            variable Z2(N,L);
        end

        Z = Z1+Z2;
        minimize( norm_nuc(Z1) + 0.1*norm_nuc(Z2) + taux*wx'*norms(Z,2,2) + ...
                  0.1*tauh*norms(Z1, 2, 1)*wh1' + tauh*norms(Z2, 2, 1)*wh2');
        
        subject to
            B*Z(:) == y;
            for kk = 1:N_known-1
                error('We should not be doing known indexes now.')
                Z(known_indexes(kk),:)*x_known(kk+1) == x_known(kk)*Z(known_indexes(kk+1),:);
            end
    cvx_end

    if isempty(strfind(cvx_status, 'Solved'))
        fname = sprintf('failed_problem_sparse_bss_nuclear_v%s', ...
                        datestr(now, 'ddmmyyyyHHMMSS'));
        warning(sprintf('cvx_status not Solved, saving %s.', fname))
        save(fname)
    end

    difference = norm(Z - Z_old,'fro')/norm(Z_old,'fro');
    
    if ~isempty(strfind(cvx_status,'Infeasible'))
        % Stop the algorithm
        if verbose
            fprintf('Infeasible cvx\n')
        end
        flag = 0;
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
