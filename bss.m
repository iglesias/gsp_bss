function bss(taux, tauh, verbose)

if nargin < 1
  taux = 1e-1;
end

if nargin < 2
  tauh = 1e-1;
end

if nargin < 3
  verbose = false;
end

bss_gen_problem

sparse_bss_verbose = false;
[Z1_hat, Z2_hat] = sparse_bss_nuclear(y, A, V, taux, tauh, sparse_bss_verbose, knownSupportFlag);
%[Z1_hat, Z2_hat] = sparse_bss_logdet(y, A, V, taux, tauh, sparse_bss_verbose, knownSupportFlag);
%[Z1_hat, Z2_hat] = sparse_bss_logdet_jointsum(y, A, V, 1e-1, 0, sparse_bss_verbose, knownSupportFlag);

%fprintf('%.4d %.4d\n', norm(Z1_hat, 'fro'), norm(Z2_hat, 'fro'))

%[Uz1, Sz1, Vz1] = svd(Z1');
%h1FromZ1 = sqrt(Sz1(1,1))*Uz1(:,1);
%x1FromZ1 = sqrt(Sz1(1,1))*Vz1(:,1);

if verbose
  bss_print_summary
%  bss_plot_results
end

end
