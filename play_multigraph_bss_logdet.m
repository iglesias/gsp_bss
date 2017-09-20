N = 100;
success = zeros(N, 1);
verbose_multigraph_bss_logdet = false;

ppm = ParforProgMon('Work', N);
parfor n = 1:N
  [truth, model, y] = multigraph_bss_gen_problem;
  Z_hat = multigraph_bss_logdet(y, model.A, model.V, ...
                                verbose_multigraph_bss_logdet);
  if recovery_assessment(truth.Z, Z_hat) < 1e-3
    success(n) = 1;
  end
  ppm.increment();
end

sum(success)/N
