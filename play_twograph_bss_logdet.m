N = 50;
success = zeros(N, 1);
verbose_twograph_bss_logdet = false;

ppm = ParforProgMon('Work', N);
parfor n = 1:N
  [truth, model, y] = twograph_bss_gen_problem;
  [Z1_hat, Z2_hat] = twograph_bss_logdet(y, model.A1, model.G(1).V, model.A2, model.G(2).V, verbose_twograph_bss_logdet);
  if recovery_assessment(truth.Z1, truth.Z2, Z1_hat, Z2_hat) < 1e-3
    success(n) = 1;
  end
  ppm.increment();
end

sum(success)/N
