% 18 September 2022
% 19

clearvars

num_mc = 50;
e = zeros(num_mc, 1);

params.S = 1;
params.L = 5;

params.filterDistribution = DataDistribution.HeatKernel;
parfor i = 1:num_mc
  [truth, model, y] = multigraph_bss_gen_problem(params);
  Z_hat = multigraph_bss_nuclear_direct(y, model.A, model.V);
  e(i) = recovery_assessment(truth.Z, Z_hat);
end
sum(e<0.1)

params.filterDistribution = DataDistribution.Uniform;
parfor i = 1:num_mc
  [truth, model, y] = multigraph_bss_gen_problem(params);
  Z_hat = multigraph_bss_nuclear_direct(y, model.A, model.V);
  e(i) = recovery_assessment(truth.Z, Z_hat);
end
sum(e<0.1)

params.filterDistribution = DataDistribution.Normal;
parfor i = 1:num_mc
  [truth, model, y] = multigraph_bss_gen_problem(params);
  Z_hat = multigraph_bss_nuclear_direct(y, model.A, model.V);
  e(i) = recovery_assessment(truth.Z, Z_hat);
end
sum(e<0.1)