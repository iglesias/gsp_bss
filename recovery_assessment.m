function v = recovery_assesment(Z1, Z2, Z1_hat, Z2_hat)

v = sqrt(norm(Z1_hat-Z1, 'fro')^2 + norm(Z2_hat-Z2, 'fro')^2) / ...
    sqrt(norm(Z1, 'fro')^2 + norm(Z2, 'fro')^2);

end
