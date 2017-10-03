function wrapper_bss_logdet_svd

verbose_bss_logdet = true;
verbose_self = true;

[truth, model, y] = bss_gen_problem;
[Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);

[UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
numFilters = length(truth.Z);
Z_hat = zeros([size(Zsum_hat) numFilters]);
for i = 1:numFilters
  Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
end

if verbose_self
  bss_svd_print_summary(Z_hat, truth, model, y);
end

end
