function multigraph_bss_print_summary(Z_hat, truth, model, y)

svd_fmt_str = prepare_svd_fmt_str(size(Z_hat, 2));

P = size(Z_hat, 3);

y_hat = 0;
for i = 1:P
  y_hat = y_hat + model.G(i).V*model.A{i}*vec(Z_hat(:, :, i));
end

fprintf('\n\n')
fprintf('Equality constraint test: %d\n', norm(y - y_hat))

for i = 1:P
  fprintf('norm(Z%d_hat-Z%d)=%d\n', i, i, norm(Z_hat(:, :, i) - truth.Z{i}))
end

for i = 1:P
  fprintf(sprintf('svd(Z%d_hat)=(%s)\n', i, svd_fmt_str), ...
          svd(Z_hat(:, :, i)))
end

fprintf('Recovery assessment: %d\n', ...
        recovery_assessment(truth.Z, Z_hat))

fprintf('\n\n')

end
