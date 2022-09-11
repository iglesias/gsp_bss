function v = recovery_assessment(Z, Z_hat)
%RECOVERY_ASSESSMENT
%   v = recovery_assessment(Z, Z_hat)
%

assert(length(Z) == size(Z_hat, 3))

% Shorthand.
R = length(Z);

v = 0;
for i = 1:R
  v = v + 1/R*norm(Z_hat(:, :, i) - Z{i}, 'fro');
end

end
