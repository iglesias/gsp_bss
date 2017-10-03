function bss_svd_print_summary(Z_hat, truth, model, y)

svd_fmt_str = prepare_svd_fmt_str(size(Z_hat, 2));

numFilters = size(Z_hat, 3);

y_hat = model.G.V*model.A*vec(sum(Z_hat, 3));

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', norm(y - y_hat))
fprintf('norm(Z-Z_hat)=%d\n', norm(truth.Zsum - sum(Z_hat, 3), 'fro') / ...
                              norm(truth.Zsum, 'fro'));

for i = 1:numFilters
  fprintf(sprintf('svd(Z%d)=(%s)\n', i, svd_fmt_str), svd(truth.Z{i}))
end

for i = 1:numFilters
  fprintf(sprintf('svd(Z%d_hat)=(%s)\n', i, svd_fmt_str), ...
          svd(Z_hat(:, :, i)))
end

fprintf('Recovery assessment: %d\n', ...
        recovery_assessment_perms(truth.Z, Z_hat));

fprintf('\n\n')

end
