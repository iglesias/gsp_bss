function wrapper_singlegraph_bss_logdet

verbose_bss_logdet = false;
verbose_self = true;

params.S = 2;
perf_sum = inf;
while perf_sum > 1e-3
  [truth, model, y] = singlegraph_svd_bss_gen_problem(params);
  [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);
  perf_sum = norm(truth.Zsum - Zsum_hat, 'fro') / norm(truth.Zsum, 'fro');
end

[UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
numFilters = length(truth.Z);
Z_hat = zeros([size(Zsum_hat) numFilters]);
for i = 1:numFilters
  Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
end

if verbose_self
  do_perms = true;
  singlegraph_bss_print_summary(Z_hat, truth, model, y, do_perms);
  plot_Zs(truth.Z, Z_hat)
  print_analysis
end

end
