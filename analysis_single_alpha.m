function data = analysis_single_alpha(alpha_arg)

dir_prefix = 'data/alpha';
fname_pattern = sprintf('%s/bss_nuclear_alpha=%.2f_10_01_*', dir_prefix, alpha_arg);

files = dir(fname_pattern);
assert(length(files) == 100)

data.eq_constraint_test = zeros(length(files), 1);
data.Zsum_test = zeros(length(files), 1);
data.Zi_norms = zeros(length(files), 2);
data.Zi_hat_norms = zeros(length(files), 2);

for i = 1:length(files)
  load(sprintf('%s/%s', dir_prefix, files(i).name), '-mat')

  data.eq_constraint_test(i) = norm(y - V*A*(Z1_hat(:) + Z2_hat(:)));
  data.Zsum_test(i) = norm([Z1+Z2] - [Z1_hat+Z2_hat], 'fro');
  data.Zi_norms(i, :) = [norm(Z1, 'fro') norm(Z2, 'fro')];
  data.Zi_hat_norms(i, :) = [norm(Z1_hat, 'fro') norm(Z2_hat, 'fro')];
end

end
