function singlegraph_bss_print_summary(Z_hat, truth, model, y, do_perms)

if ~exist('do_perms', 'var')
  do_perms = false;
end

svd_fmt_str = prepare_svd_fmt_str(size(Z_hat, 2));

P = size(Z_hat, 3);

y_hat = model.G.V*model.A*vec(sum(Z_hat, 3));

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', norm(y - y_hat))
fprintf('norm(Z-Z_hat)=%d\n', norm(truth.Zsum - sum(Z_hat, 3), 'fro') / ...
                              norm(truth.Zsum, 'fro'));

for i = 1:P
  fprintf('norm(Z%d_hat-Z%d)=%d\n', i, i, norm(Z_hat(:, :, i) - truth.Z{i}))
end

for i = 1:P
  fprintf(sprintf('svd(Z%d)=(%s)\n', i, svd_fmt_str), svd(truth.Z{i}))
end

for i = 1:P
  fprintf(sprintf('svd(Z%d_hat)=(%s)\n', i, svd_fmt_str), ...
          svd(Z_hat(:, :, i)))
end

if do_perms
  perf = recovery_assessment_perms(truth.Z, Z_hat);
else
  perf = recovery_assessment(truth.Z, Z_hat);
end

fprintf('Recovery assessment: %d\n', perf);

fprintf('\n\n')

end
