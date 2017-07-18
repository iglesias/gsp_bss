clearvars

alphas = 0:0.1:1;
num_mc = 100;

for alpha = alphas
  parfor n = 1:num_mc
    bss(alpha)
  end

  fprintf('Alpha %.2f done!\n', alpha)
end
