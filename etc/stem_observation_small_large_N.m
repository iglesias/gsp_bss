clearvars

params.L = 2;
params.N = 100;
params.S = 1;
params.numGraphs = 5;

max_small_N = 10;
max_large_N = 10;

while true
  params.N = 30;
  [truth, model, y_small_N] = multigraph_bss_gen_problem(params);
  max_small_N = max(max_small_N, max(y_small_N));
  subplot(211)
  stem(y_small_N)
  axis([1 params.N 0 max_small_N])
  
  params.N = 100;
  [truth, model, y_large_N] = multigraph_bss_gen_problem(params);
  max_large_N = max(max_large_N, max(y_large_N));
  subplot(212)
  stem(y_large_N)
  axis([1 params.N 0 max_large_N])
  
  pause(0.5)
end