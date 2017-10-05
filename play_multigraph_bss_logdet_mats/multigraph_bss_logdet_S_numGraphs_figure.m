% Edit the last element to select number of columns / values of S.
SS = 1:8;
NUM_GRAPHS = 2:5;

data_matrix = [0.97 0.92 0.77 0.67; ...
               0.91 0.80 0.49 0.30; ...
               0.89 0.61 0.34 0.08; ...
               0.71 0.31 0.10 0.00; ...
               0.64 0.26 0.02 0.00; ...
               0.52 0.20 0.00 0.00; ...
               0.58 0.07 0.00 0.00; ...
               0.38 0.03 0.00 0.00]';

figure
imagesc(data_matrix)
axis xy
xticks(SS)
yticks(NUM_GRAPHS-1)
yticklabels(NUM_GRAPHS)
xlabel('S')
ylabel('Number of graphs')
colormap gray
colorbar
title('Rate of successful recovery with L=3 and N=100')
