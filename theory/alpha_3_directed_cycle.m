clearvars

normalize_Psi = true;

% Number of nodes.
% NN = [100 200 400 800 1600];
NN = [100 300 600];
alpha_3 = zeros(length(NN), 1);
for n = 1:length(NN)
  N = NN(n);
  % Adjacency matrix.
  model.G.W = circshift(eye(N), [N 1]); % directed cycle
%   p = 0.5;
%   model.G.W = generate_connected_ER(N, p);

  % Filter coefficients.
  L = 1;
  assert(L <= N)
  h = randn(L, 1);
  % Normalize filter taps.
  h = h / norm(h);

  [model.G.V, Lambda] = eig(model.G.W);
  model.G.U = inv(model.G.V);
  U = model.G.U;

  model.G.lambda = diag(Lambda);
  model.Psi = repmat(model.G.lambda, 1, L).^repmat([0:L-1], N, 1); %#ok<NBRAK>
  assert(N == size(model.Psi, 1))

  if normalize_Psi
    [Psi_svd_L, ~, ~] = svd(model.Psi, 0);
    Psi = Psi_svd_L;
  else
    Psi = model.Psi;
  end

  mu_max_squared = max(N/L*max(norms(Psi, 2, 2).^2));
%   display(mu_max_squared)
  mu_h_squared = N*max(max(abs(conj(Psi)*h).^2) / norm(h)^2);
%   display(mu_h_squared)

  R = 2;

  % % % Gaussian model demixing
%   alpha_3(n) = N / (R * max(mu_max_squared*L, mu_h_squared*N) * (log(N))^2) - log(R);

  % % % GSP demixing.
  % First submission:
%   alpha_3(n) = 1 / (N*R * max(mu_max_squared*L, mu_h_squared*N^2) * (log(N))^2) - log(R); 
  % After review:
  S = 2;
  assert(S <= N)
  rho_max = rhof(U, 1);
  alpha_3(n) = N / (R*(rho_max + 1)^2 * max(mu_max_squared*L, mu_h_squared*S^2) * (log(N))^2) - log(R); 
  
%   rho_U__1 = max(max(abs(U).^2));
%   max_second_arg = mu_h_squared*N*(sqrt(rho_U__1*N) + 1)^2;
%   alpha_3(n) = N / (R * max(mu_max_squared*L, max_second_arg) * (log(N))^2) - log(R);
end

figure
plot(NN, alpha_3, 'o-', 'LineWidth', 1.5)
grid on
xlabel('N')
ylabel('\alpha_3')

if normalize_Psi
  title('Vandermonde normalized')
else
  title('Vandermonde w/o normalization')
end