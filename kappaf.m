function kappa = kappaf(A, B, k1, k2)

assert(all(size(A) == size(B)))
assert(size(A, 2) >= k1)
assert(size(B, 2) >= k2)

kappa = -Inf;
for ell = 1:size(A, 1)
  % Non-vectorized computation of
  % max_{(Omega_1, Omega_2) \in ...} ||a_{ell_{Omega_1}} b_{ell_{Omega_2}}^H||
  X = nchoosek(A(ell, :), k1);
  Y = nchoosek(B(ell, :), k2);
  for i = 1:size(X, 1)
    for j = 1:size(Y, 1)
      kappa = max(kappa, norm( X(i,:).' * conj(Y(j,:)) ));
    end
  end 
end

end