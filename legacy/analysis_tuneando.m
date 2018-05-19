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

    pause(5)
end

end
