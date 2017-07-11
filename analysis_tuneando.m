function analysis_tuneando

fname_prefix = 'data/tuneando/bss_nuclear_tuneando_10_01*';

files = dir(fname_prefix);

for i = 1:length(files)
    load(sprintf('data/tuneando/%s', files(i).name))

    cminmax = minmax([Z1(:); Z2(:); Z1_hat(:); Z2_hat(:)]);

    subplot(221)
    imagesc(Z1, cminmax), title('Z1')
    subplot(222)
    imagesc(Z2, cminmax), title('Z2')
    subplot(223)
    imagesc(Z1_hat, cminmax), title('Z1\_hat')
    subplot(224)
    imagesc(Z2_hat, cminmax), title('Z2\_hat')

    clc
    fprintf('%s (%d/%d)\n', files(i).name, i, length(files))
    fprintf('\n\n')
    fprintf('Equality constraint test: %d\n', norm(y - V*A*(Z1_hat(:) + Z2_hat(:))))
    fprintf('norm(Z1)=%d\n', norm(Z1))
    fprintf('norm(Z2)=%d\n', norm(Z2))
    fprintf('norm(Z1_hat)=%d\n', norm(Z1_hat))
    fprintf('norm(Z2_hat)=%d\n', norm(Z2_hat))
    fprintf('norm(Z1_hat-Z1)=%d\n', norm(Z1_hat - Z1))
    fprintf('norm(Z2_hat-Z2)=%d\n', norm(Z2_hat - Z2))
    fprintf('norm(Z1_hat-Z2)=%d\n', norm(Z1_hat - Z2))
    fprintf('norm(Z2_hat-Z1)=%d\n', norm(Z2_hat - Z1))
    fprintf('norm([Z1+Z2]-[Z1_hat+Z2_hat])=%d\n', norm([Z1+Z2] - [Z1_hat+Z2_hat]))
    fprintf('norm(Z1-Z2)=%d\n', norm(Z1 - Z2))
    fprintf('norm(Z1_hat-Z2_hat)=%d\n', norm(Z1_hat - Z2_hat))
    fprintf('norm_nuc(Z1)+norm_nuc(Z2)=%d\n', norm_nuc(Z1)+norm_nuc(Z2))
    fprintf('norm_nuc(Z1_hat)+norm_nuc(Z2_hat)=%d\n', norm_nuc(Z1_hat)+norm_nuc(Z2_hat))
    pause(5)
end

end
