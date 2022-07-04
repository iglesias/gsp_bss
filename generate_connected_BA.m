function A = generate_connected_BA(N, m)

assert(N > m)

is_connected = false;
B = generate_connected_ER(m, 0.5);

while ~is_connected

A = zeros(N);
A(1:m, 1:m) = B;

for i = (m+1):N
    d = sum(A);
    if ~all(d((i+1):end) == 0)
        keyboard
    end
    l = 0;
    while l < min(m, sum(find(A(i, 1:(i-1)) == 0)))
        j = find(cumsum(d)/sum(d) > rand, 1);
        if isempty(j), assert(false), end
        assert(j <= (i-1))
        if A(i,j) == 0
            A(i,j) = 1;
            A(j,i) = 1;
            l = l+1;
        end
    end
end

% Is connected?
L = diag(sum(A)) - A;
[~, D] = eig(L);
if sum(abs(diag(D)) < 1e-5) == 1
    is_connected = true;
end

end

end
