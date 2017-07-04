function [Z1_hat, Z2_hat] = sparse_bss_nuclear(y,A,V,taux,tauh,varargin)
% SPARSE_BSS_NUCLEAR: Obtains the estimate of the filter and the elements 
% of a graph signal using the nuclear norm surrogate, assuming a sparse model for the signal
%
%
%           [Z1_hat, Z2_hat] = sparse_bss_nuclear(y,M,L,taux,tauh,x_known,known_indexes)
%           y = filtered graph signal (might be a subset)
%           A = Khatri-Rao product matrix
%           V = rows of the inverse Graph Fourier transform which correspond to the observations
%           taux = regularization parameter for x
%           tauh = regularization parameter for h
%           knownSupportFlag = whether the support of the input signals (assumed the same) is known
%           S = cardinality of the input signals' support
%           x_known = known values of the graph signal x. Could be an empty matrix or optional
%           known_indexes = indexes of the known values. Could be an empty matrix or optional
%       
%           Zi_hat = estimate of the rank-1 matrices xi*hi'

% David Ramirez
% Signal Processing group
% University Carlos III of Madrid
% 2016

%% Disable Warnings

warning ('off','MATLAB:nargchk:deprecated')

%% Parameter definition

epsilon_normx = 1e-3;
epsilon_normh = 1e-3;
maxiter = 50;

%% Check input parameters

if nargin<7
    knownSupportFlag = 0;
    S = [];
else 
    assert(nargin==7 || nargin==9)
    knownSupportFlag = varargin{1};
    S = varargin{2};
end

if nargin<9
    known_indexes = [];
    x_known = [];
else
    assert(nargin==9)
    x_known = varargin{3};
    known_indexes = varargin{4};
end

%% Blind estimation using rank minimization: nuclear norm surrogate

assert(size(A, 1) == length(y)) % this might change for sampled observations
N = size(A, 1);
% It is assumed that the filters have the same order.
if knownSupportFlag
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
if knownSupportFlag
    Z_old = reshape(pinv(B)*y,S,L);
else
    Z_old = reshape(pinv(B)*y,N,L);
end

while (flag == 1 && iter <= maxiter)
    fprintf('iteration %d\n', iter)

    wx = 1./(sqrt(sum(abs(Z_old).^2,2)) + epsilon_normx);
    wx(known_indexes) = 0;
    wh = 1./(sqrt(sum(abs(Z_old).^2,1)) + epsilon_normh);
    
    cvx_begin quiet
        if knownSupportFlag
            variable Z1(S,L);
            variable Z2(S,L);
        else
            variable Z1(N,L);
            variable Z2(N,L);
        end

        Z = Z1+Z2;
        minimize( norm_nuc(Z1) + norm_nuc(Z2) + taux*wx'*norms(Z,2,2) + tauh*norms(Z,2,1)*wh' );
        
        subject to
            B*Z(:) == y;
            for kk = 1:N_known-1
                Z(known_indexes(kk),:)*x_known(kk+1) == x_known(kk)*Z(known_indexes(kk+1),:);
            end
    cvx_end

    if ~strcmp(cvx_status, 'Solved')
        warning('cvx_status not equal to Solved.')
        keyboard
    end

    difference = norm(Z - Z_old,'fro')/norm(Z_old,'fro');
    
    if ~isempty(strfind(cvx_status,'Infeasible'))
        % Stop the algorithm
        fprintf('Infeasible cvx\n')
        Z = randn(N,L);
        flag = 0;
    else
        if difference<1e-4
            % Converged
            fprintf('Convergence reached\n')
            flag = 0;
        else
            % Did not converge
            fprintf('Covergence NOT reached, difference=%d.\n', difference)
            Z_old = Z;
            iter = iter + 1;
        end
    end
        
end

Z1_hat = Z1;
Z2_hat = Z2;
