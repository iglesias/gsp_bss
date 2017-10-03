function v = recovery_assessment_perms(Z, Z_hat)

assert(length(Z) == size(Z_hat, 3))

numFilters = size(Z_hat, 3);
v = Inf;
P = perms([1:numFilters]);
for p = 1:size(P, 1)
  for i = 1:numFilters
    Z_hat_rearranged(:, :, i) = Z_hat(:, :, P(p, i));
  end

  v = min(recovery_assessment(Z, Z_hat_rearranged), v);
end

end
