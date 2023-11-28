% Modified Dec 21 2020

clearvars
N = 102;
num_mc = 5e2;

%% |<u_{1,l}, x_1>|
R = zeros(num_mc, N);

for i = 1:num_mc
    x = randn(N, 1);
    x = x / norm(x);
    assert(abs(norm(x) - 1) < 1e-5)

    S = generate_connected_ER(N, 0.1);
    [V, ~] = eig(S);
    U = inv(V);
    U = sqrt(N)*U;

    assert(norm(U'*U - N*eye(N)) < 1e-5)
    
    R(i, :) = U*x;
end

figure
histogram(R(:))
title(sprintf('Mean of |samples|^2: %.5f', mean(R(:).^2)))

%% ||conj(u_{2,l}) x_1'||

T = zeros(num_mc, 1);

for i = 1:num_mc
    x = randn(N, 1);
    x = x / norm(x);
    assert(abs(norm(x) - 1) < 1e-5)

    S = generate_connected_ER(N, 0.1);
    [V, ~] = eig(S);
    U = inv(V);
    U = sqrt(N)*U;
    assert(norm(U'*U - N*eye(N)) < 1e-5)
    
    assert(isreal(x))
    T(i) = norm(conj(U(2,:)'*x'));    
end

figure
histogram(T(:))
title(sprintf('sqrt(N)=%.5f min=%.5f max=%.5f', sqrt(N), min(T(:)'), max(T(:)')))

%%
clc
L = 10;
max_value = 0;
A = 0;
for i = 1:5
    v = randn(L, 1) + 1i*randn(L, 1);
    A = A + conj(v)*v';
    max_value = max(max_value, norm(v)^2);
end

max_value
norm(A)
svd(A)
