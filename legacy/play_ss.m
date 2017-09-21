clearvars

% Number of nodes.
N = 10;

% Adjacency matrix.
G = generate_connected_ER(N, 0.1); 

% Filters.
num_filters = 2;

L1 = [0.9 0.1 0];
L2 = [0.95 0 0.05];

H1 = L1(1)*eye(N) + L1(2)*G + L1(3)*G*G;
H2 = L2(1)*eye(N) + L2(2)*G + L2(3)*G*G;

H = [H1 H2];

% Inputs.
x1 = zeros(N, 1);
x2 = zeros(N, 1);

x1(1:2) = randn(2, 1);
x2(1:2) = randn(2, 1);

x_true = [x1; x2];

% Output.
y = H*x_true;

% Source separation.

% Least-squares solution to initialize the estimate.
xhat_old = pinv(H)*y;

epsilon = 1e-5;
iter_idx = 1;
max_iter = 10;
continue_flag = true;

while iter_idx <= max_iter && continue_flag 
  fprintf('Iteration %d...\n', iter_idx)
  w = 1 ./ (abs(xhat_old) + epsilon);
  cvx_begin quiet
    variable x(num_filters*N);

    minimize(w'*abs(x));
    subject to
      y == H*x;
  cvx_end

  xhat = x;
  difference = norm(xhat - xhat_old) / norm(xhat_old);
  fprintf('  difference=%d\n', difference)

  if difference < 1e-6
    fprintf('Converged after %d iterations\n', iter_idx);
    continue_flag = false;
  end

  xhat_old = xhat;
  iter_idx = iter_idx + 1;
end
