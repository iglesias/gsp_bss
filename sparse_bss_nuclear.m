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

if nargin==5
    known_indexes = [];
    x_known = [];
else
    x_known = varargin{1};
    known_indexes = varargin{2};
end

%% Blind estimation using rank minimization: nuclear norm surrogate

assert(size(A, 1) == length(y)) % this might change for sampled observations
N = size(A, 1);
% It is assumed that the filters have the same order.
L = size(A, 2)/N;
B = V*A;

% Length of the known signal values
N_known = length(known_indexes);

% Initialization
flag = 1;
iter = 1;
Z_old = reshape(pinv(B)*y,N,L);         % LS solution for initializing w

while (flag == 1 && iter <= maxiter)

    wx = 1./(sqrt(sum(abs(Z_old).^2,2)) + epsilon_normx);
    wx(known_indexes) = 0;
    wh = 1./(sqrt(sum(abs(Z_old).^2,1)) + epsilon_normh);
    
    cvx_begin quiet
        variable Z1(N,L);
        variable Z2(N,L);

        Z = Z1+Z2;
        minimize( norm_nuc(Z) + taux*wx'*norms(Z,2,2) + tauh*norms(Z,2,1)*wh' );
        
        subject to
            B*Z(:) == y;
            for kk = 1:N_known-1
                Z(known_indexes(kk),:)*x_known(kk+1) == x_known(kk)*Z(known_indexes(kk+1),:);
            end
    cvx_end

    difference = norm(Z - Z_old,'fro')/norm(Z_old,'fro');
    
    if ~isempty(strfind(cvx_status,'Infeasible'))
        % Stop the algorithm
        Z = randn(N,L);
        flag = 0;
    else
        if difference<1e-4
            % Converged
            flag = 0;
        else
            % Did not converge
            Z_old = Z;
            iter = iter + 1;
        end
    end
        
end

Z1_hat = Z1;
Z2_hat = Z2;
