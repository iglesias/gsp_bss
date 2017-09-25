function test_input_variation

N = 100;
svd_mean = 0;
l11_norm_mean = 0;
l21_norm_mean = 0;

for i = 1:N
  [truth, model, y] = multigraph_bss_gen_problem;
  svd_mean = svd_mean + svd(truth.Z{1})/N;
  l11_norm_mean = l11_norm_mean + sum(norms(truth.Z{1}, 1, 2))/100;
  l21_norm_mean = l21_norm_mean + sum(norms(truth.Z{1}, 2, 2))/100;
end

display(svd_mean)
display(l11_norm_mean)
display(l21_norm_mean)

end
