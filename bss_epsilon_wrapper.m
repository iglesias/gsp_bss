clearvars

epsilons = [1e-10 1e-8];
num_mc = 100;

for epsilon = epsilons
  parfor n = 1:num_mc
    bss(epsilon, 0.1, [], false)
  end

  fprintf('Epsilon %.0d done!\n', epsilon)
end
