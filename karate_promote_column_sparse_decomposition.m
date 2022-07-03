N = 34;

x1 = zeros(N, 1);
x2 = zeros(N, 1);

todo = true;
while todo
  x1_support = randperm(N, 1);
  x2_support = randperm(N, 2);
  
  if x1_support ~= x2_support
    todo = false;
  end
end

x1(x1_support) = 1;
x2(x2_support) = 1;

assert(abs(x1'*x2) < 1e-7)

% orthogonal 2d-vectors of filter taps.
% h1' = [h11 h12], h2' = [h21 h22]
% h1'*h2 = 0 => h11*h21 + h12*h22 = 0
h1 = randn(2, 1);
h21 = randn;
h22 = -h1(1)*h21/h1(2);
h2 = [h21 h22]';

h1 = h1 / norm(h1);
h2 = h2 / norm(h2);

assert(abs(h1'*h2) < 1e-7)
L = 2;

Z1 = x1*h1';
Z2 = x2*h2';
Zsum = Z1 + Z2;

cvx_solver mosek
cvx_begin quiet
  variable Z1_hat(N, L);
  variable Z2_hat(N, L);
  
  minimize(norm(sum(Z1_hat, 2), 1) + norm(sum(Z2_hat, 2), 1) + norm_nuc(Z1_hat) + norm_nuc(Z2_hat))
  
  subject to
    Z1_hat + Z2_hat == Zsum;
cvx_end

Z_hat = zeros(N, L, 2);
Z_hat(:, :, 1) = Z1_hat;
Z_hat(:, :, 2) = Z2_hat;
do_perms = true;
truth.Z = cell(2, 1);
truth.Z{1} = Z1;
truth.Z{2} = Z2;
truth.Zsum = Zsum;
singlegraph_bss_print_summary(Z_hat, truth, [], [], do_perms, Z1_hat+Z2_hat);
plot_Zs(truth.Z, Z_hat)