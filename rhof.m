function rho = rhof(A, k)

assert(size(A, 2) >= k)

rho = -Inf;
for ell = 1:size(A, 1)
  rho = max(rho, max(norms(nchoosek(A(ell, :), k), 2, 2).^2));
end

end

