clear vars

v = [1 2]';
w = [4 5]';
X = kron(v, w);
Xtilda = [];
i = 1; j = 2;
while j <= length(X)
  Xtilda = [Xtilda; X(i:j)'];
  i = j+1;
  j = j+2;
end

[U,S,V] = svd(Xtilda);
