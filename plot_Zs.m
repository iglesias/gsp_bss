function plot_Zs(Z, Z_hat)
  
if length(Z) ~= size(Z_hat, 3)
  error('Size mismatch between Z and Z_hat')
end

num_subplot_cols = length(Z);
for i = 1:num_subplot_cols
  subplot(2, num_subplot_cols, i)
  imagesc(Z{i})
  subplot(2, num_subplot_cols, i + num_subplot_cols)
  imagesc(Z_hat(:, :, i))
end

end
