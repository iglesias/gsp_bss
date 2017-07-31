cminmax = minmax([Z1(:); Z2(:); Z1_hat(:); Z2_hat(:)]');

figure(1)
clf
subplot(221)
imagesc(Z1, cminmax), title('Z1')
subplot(222)
imagesc(Z2, cminmax), title('Z2')
subplot(223)
imagesc(Z1_hat, cminmax), title('Z1\_hat')
subplot(224)
imagesc(Z2_hat, cminmax), title('Z2\_hat')

clear cminmax
