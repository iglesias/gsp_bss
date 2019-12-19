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
indices = find(NN == 20);
i = indices(end); j = indices(end-1);

y = XX{i} + XX{j};

Psi_cut = 10;

tmp_Psi_i = fliplr(vander(EE{i}));
Psi_i = tmp_Psi_i(:, 1:Psi_cut);
U_i = VV{i}';
model.A{1} = kr(Psi_i', U_i')';
tmp_Psi_j = fliplr(vander(EE{j}));
Psi_j = tmp_Psi_j(:, 1:Psi_cut);
U_j = VV{j}';
model.A{2} = kr(Psi_j', U_j')';

model.V = zeros(20, 20, 2); % two graphs/sources of 20 nodes.
model.V(:,:,1) = VV{i};
model.V(:,:,2) = VV{j};

verbose_solver = true;
Z_hat = multigraph_bss_logdet(y, model.A, model.V, verbose_solver);

model.G(1).V = model.V(:,:,1);
model.G(2).V = model.V(:,:,2);
truth.Z = {SS{i}*(pinv(Psi_i)*HH_Freq{i})', SS{j}*(pinv(Psi_j)*HH_Freq{j})'};
multigraph_bss_print_summary(Z_hat, truth, model, y);
