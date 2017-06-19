function [x_hat] = knownh_bci(y,H)
% KNOWNH_BCI: Obtains the estimate of the elements of a graph 
% signal assuming the filter known
%
%
%           [x_hat] = knownh_bci(y,H)
%           y = filtered graph signal (might be a subset)
%           H = rows of the graph filter that correspond to the observations
%       
%           x_hat = estimate of the graph signal

% David Ramirez
% Signal Processing group
% University Carlos III of Madrid
% 2017

%% Disable Warnings

warning ('off','MATLAB:nargchk:deprecated')

%% Linopt options

options = optimoptions(@linprog,'Algorithm','interior-point','MaxIter',200,'TolFun',1e-8,'Display','off','Diagnostics','off');

%% Parameter definition

epsilon = 1e-3;
maxiter = 50;

%% Blind estimation using rank minimization: nuclear norm surrogate

% Obtain the number of nodes
[~,N] = size(H);

% Initialization
flag = 1;
iter = 1;
x_old = pinv(H)*y;         % LS solution for initializing w

while (flag == 1 && iter <= maxiter)

    w = 1./(abs(x_old) + epsilon);
    
    % cvx_begin quiet
    %     variable x(N);
    %     
    %     minimize( w'*abs(x) )
    %     subject to
    %         H*x == y;
    % cvx_end
    
        [z,~,exitflag] = linprog([w ; w],[],[],[H -H],y,zeros(2*N,1),[],[],options);
    
    if (exitflag<=0)
        % Stop the algorithm
        x_hat = randn(N,1);
        flag = 0;
    else
        x_hat = z(1:N) - z(N+1:end);
        
        difference = norm(x_hat - x_old)/norm(x_old);
        
        if difference<1e-4
            % Converged
            flag = 0;
        else
            % Did not converge
            x_old = x_hat;
            iter = iter + 1;
        end
    end
        
end
