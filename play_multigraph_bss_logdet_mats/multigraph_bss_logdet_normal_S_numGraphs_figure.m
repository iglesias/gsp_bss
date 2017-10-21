% Edit the last element to select number of columns / values of S.
SS = 1:5;
NUM_GRAPHS = 2:5;

data_matrix = [0.91 0.79 0.75 0.57; ...
               0.84 0.54 0.27 0.15; ...
               0.74 0.33 0.12 0.02; ...
               0.65 0.26 0.03 0.00; ...
               0.51 0.15 0.01 0.00]';

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
title('Rate of successful recovery with L=3 and N=100')
