clearvars

% Load data from dataset
load('dataset\all_data_eigL_x6');

% Control variables
n_graphs = min(100, length(XX));
NN = zeros(n_graphs, 1);
max_gso_res = zeros(n_graphs, 1);

GSO = cell(n_graphs, 1);

for i = 1:n_graphs
    % Select data from graph i
    h_freq = HH_Freq{i};
    e = EE{i};
    V = VV{i};
    x = XX{i};
    s = SS{i};
    NN(i) = length(h_freq);
    Psi = fliplr(vander(e));
    H = V*diag(h_freq)*V';   
%     disp(['Graph: ' num2str(i) '   N_nodes: ' num2str(NN(i))])
    
    GSO{i} = V*diag(e)*V';
    max_gso_res(i) = max(vec(abs(GSO{i} - round(GSO{i}))));
end

% There are 10 graphs with 20 nodes.
