% Edit the last element to select number of columns / values of S.
SS = 1:6;
NUM_GRAPHS = 2:6;

data_matrix = [0.829 0.668 0.479 0.283 0.126; ...
               0.603 0.342 0.110 0.022 0.000; ...
               0.455 0.145 0.020 0.001 0.000; ...
               0.301 0.043 0.001 0.000 0.000; ...
               0.198 0.009 0.000 0.000 0.000]';

figure
imagesc(data_matrix)
axis xy
xticks(SS)
yticks(NUM_GRAPHS-1)
yticklabels(NUM_GRAPHS)
xlabel('S')
ylabel('Number of graphs')
caxis([0 1])
colormap gray
colorbar
title('Rate of successful recovery with brain graphs (N=66, L=3)')
