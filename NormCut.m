function [z stats] = NormCut(subject, experiment, n_clust)

addpath('normcut');
stats = struct('NMI',[], 'conn_diff',[]);
loaded = load(['../data/' subject '/' experiment '.mat']);
D = loaded.D;
adj_list = loaded.adj_list;
coords = loaded.coords;
N = size(D,1);

if (strcmp(experiment, 'PPA'))
    labels = loaded.labels;
    bold = loaded.bold;
end
if (strcmp(subject, 'synth'))
    gt_z = loaded.z;
else
    gt_z = [];
end
clear loaded;

W = zeros(size(D));
eps = 0.0001;
for i = 1:size(D,1)
    for j = adj_list{i}
        W(i,j) = 1 / (norm(D(i, (1:N ~= i) & (1:N ~= j)) - ...
                           D(j, (1:N ~= i) & (1:N ~= j)), 2) + eps);
    end
end

if (~strcmp(experiment, 'PPA'))
    nc_discrete = ncutW(W, n_clust);
    z = nc_discrete * (1:n_clust)';
    if (~isempty(gt_z))
        stats.NMI = CalcNMI(gt_z, z);
    end
else
    z = zeros(size(D,1),1);
    for side=1:2
        sideW = W(labels==(side+8),labels==(side+8));
        nc_discrete = ncutW(sideW, n_clust);
        sidez = nc_discrete*(1:n_clust)';
        z(labels==(side+8)) = sidez+n_clust*(side-1);
    end
    
    stats.conn_diff = CalcPPAConnDiff(z, labels, coords, bold);
end
end

