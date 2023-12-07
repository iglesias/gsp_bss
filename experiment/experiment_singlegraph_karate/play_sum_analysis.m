function [perf, perf_sum] = play_sum_analysis

verbose_bss_logdet = false;

num_simulations = 100;

params.L = 3;
params.numFilters = 2;
params.S = 1;

perf = zeros(num_simulations, 1);
perf_sum = zeros(num_simulations, 1);

parfor i_mc = 1:num_simulations
    [truth, model, y] = karate_bss_gen_problem(params);
    [Zsum_hat, iter] = bss_logdet_jointsum(y, model.A, model.G.V, verbose_bss_logdet);

    [UZ, SZ, VZ] = svd(Zsum_hat, 'econ');
    numFilters = length(truth.Z);
    Z_hat = zeros([size(Zsum_hat) numFilters]);
    for i = 1:numFilters
      Z_hat(:, :, i) = SZ(i,i)*UZ(:,i)*VZ(:,i)';
    end
    
    perf(i_mc) = recovery_assessment_perms(truth.Z, Z_hat);
    perf_sum(i_mc) = norm(truth.Zsum - sum(Z_hat, 3), 'fro') / norm(truth.Zsum, 'fro');
end

end
