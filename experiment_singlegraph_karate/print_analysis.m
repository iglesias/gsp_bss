Z_hat_sum = zeros(size(truth.Zsum));
Z_truth_sum = zeros(size(truth.Zsum));

for i = 1:numFilters
    Z_hat_sum = Z_hat_sum + Z_hat(:,:,i);
    Z_truth_sum = Z_truth_sum + truth.Z{i};
end

fprintf('\nsum equivalences:\n')
fprintf('%d\n', norm(truth.Zsum - Z_hat_sum, 'fro') / norm(truth.Zsum, 'fro'))
fprintf('%d\n', norm(truth.Zsum - Z_truth_sum, 'fro') / norm(truth.Zsum, 'fro'))



[UZsum, SZsum, VZsum] = svd(truth.Zsum, 'econ');
Z_svd_sum = zeros(size(truth.Zsum));
for i = 1:numFilters
    Z_svd_sum = Z_svd_sum + SZsum(i,i)*UZsum(:,i)*VZsum(:,i)';
end
fprintf('\nsvd of truth sum:\n')
fprintf('%d\n', norm(truth.Zsum - Z_svd_sum, 'fro') / norm(truth.Zsum, 'fro'))



fprintf('\nthe svd of the truth sum does not give the truth mixing components:\n')
for i = 1:numFilters
    fprintf('%d\n', norm(truth.Z{i} - SZsum(i,i)*UZsum(:,i)*VZsum(:,i)', 'fro') / norm(truth.Z{i}, 'fro'))
end



fprintf('\nranks:\n')
for i = 1:numFilters
    fprintf('%d %d ', rank(truth.Z{i}), rank(Z_hat(:,:,i)))
end
fprintf('\n')
