function v = recovery_assesment(Z, Z_hat)
%RECOVERY_ASSESMENT
%   v = recovery_assesment(Z, Z_hat)
%

assert(length(Z) == size(Z_hat, 3))

num = 0;
den = 0;
for i = 1:length(Z)
  num = num + norm(Z_hat(:, :, i) - Z{i}, 'fro')^2;
  den = den + norm(Z{i}, 'fro')^2;
end

v = sqrt(num) / sqrt(den);

end
