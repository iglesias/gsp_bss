function [x_hat,h_hat,sv] = sparse_bci_logdet(y,T,V,taux,tauh,varargin)
% SPARSE_BCI_LOGDET: Obtains the estimate of the filter and the elements of 
% a graph signal using the log-det surrogate, assuming a sparse model for the signal
%
%
%           [x_hat,h_hat,sv] = sparse_bci_logdet(y_freq,T,L,taux,tauh,x_known,known_indexes)
%           y = filtered graph signal (might be a subset)
%           T = Tensor
%           V = rows of the inverse Graph Fourier transform which correspond to the observations
%           taux = regularization parameter for x
%           tauh = regularization parameter for h
%           x_known = known values of the graph signal x. Could be an empty matrix or optional
%           known_indexes = indexes of the known values. Could be an empty matrix or optional
%       
%           x_hat = estimate of the graph signal
%           h_hat = filter coefficients
%           sv = singular values

% David Ramirez
% Signal Processing group
% University Carlos III of Madrid
% 2016

%% Disable Warnings

warning ('off','MATLAB:nargchk:deprecated')

%% Parameter definition

epsilon_rank = 5e-2;
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

%% Blind estimation using rank minimization: Log-det surrogate

% Obtain the matrix from the Tensor
[~,L,N] = size(T);
B = V*reshape(T,[L*N N])';

% Length of the known signal values
N_known = length(known_indexes);

% Initialization
flag = 1;
iter = 1;
Z1_old = eye(N);
Z2_old = eye(L);
Z_old = reshape(pinv(B)*y,N,L);         % LS solution for initializing w

% Majorization-minimization
while (flag == 1 && iter <= maxiter)

    wx = 1./(sqrt(sum(abs(Z_old).^2,2)) + epsilon_normx);
    wx(known_indexes) = 0;
    wh = 1./(sqrt(sum(abs(Z_old).^2,1)) + epsilon_normh);
    
    cvx_begin quiet
        variable Z(N,L);
        variable Z1(N,N) symmetric;
        variable Z2(L,L) symmetric;
        
        minimize( trace((Z1_old + epsilon_rank*eye(N))\Z1) + trace((Z2_old + epsilon_rank*eye(L))\Z2) + taux*wx'*norms(Z,2,2) + tauh*norms(Z,2,1)*wh' );
        subject to
            [Z1 Z ; Z' Z2] == semidefinite(N+L);
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
            Z1_old = Z1;
            Z2_old = Z2;
            Z_old = Z;
            iter = iter + 1;
        end
    end
        
end

%% Best rank-one approximation

unknown_indexes = 1:N;unknown_indexes(known_indexes) = [];

Z_K = Z(known_indexes,:);
Z_Kc = Z(unknown_indexes,:);

[V,sv] = eig(Z_Kc'*Z_Kc + Z_K'*(x_known*x_known')*Z_K,'vector');
[sv,index] = sort(sqrt(abs(sv)),'descend');
h_tilde = V(:,index(1));

if N_known>0
    scale_h_tilde = h_tilde'*Z_K'*x_known/(sum(abs(x_known.^2)));
else
    scale_h_tilde = 1;
end

h_hat = scale_h_tilde*h_tilde;
x_hat = zeros(N,1);
x_hat(known_indexes) = x_known;
x_hat(unknown_indexes) = Z_Kc*h_hat/(h_hat'*h_hat);