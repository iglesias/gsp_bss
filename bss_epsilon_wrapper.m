clearvars

epsilons = 0:0.1:1;
num_mc = 100;

for epsilon = epsilons
  parfor n = 1:num_mc
    bss(epsilon, 0.1, [], false)
  end

  fprintf('Epsilon %.2f done!\n', epsilon)
end
