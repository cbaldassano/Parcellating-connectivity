function [z stats] = LocalSimilarity(subject, experiment, n_clust)

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

W = inf(size(D));
for i = 1:size(D,1)
    for j = adj_list{i}
        W(i,j) = norm(D(i, (1:N ~= i) & (1:N ~= j)) - ...
                      D(j, (1:N ~= i) & (1:N ~= j)), 2);
    end
end

if (~strcmp(experiment, 'PPA'))
    z = ThresholdSimilarity(W, n_clust);
    stats.NMI = CalcNMI(gt_z, z);
else
    z = zeros(size(D,1),1);
    for side=1:2
        sideW = W(labels==(side+8),labels==(side+8));
        sidez = ThresholdSimilarity(sideW, n_clust);
        z(labels==(side+8)) = sidez + n_clust*(side-1);
    end
    
    stats.conn_diff = CalcPPAConnDiff(z, labels, coords, bold);
end
end

function z = ThresholdSimilarity(W, n_clust)

thresh = mean(W(isfinite(W)));
thresh_step = 0.1;
last_increase = 0;
last_decrease = 0;
while(1)
    G = sparse(W <= thresh);
    [~, z] = graphconncomp(G);
    curr_n_clust = length(unique(z));
    if (curr_n_clust == n_clust)
        break;
    elseif (curr_n_clust < n_clust)
        if (last_increase)
            thresh_step = thresh_step^2;
        end
        thresh = thresh * (1-thresh_step);
        last_decrease = 1; last_increase = 0;
    else
        if (last_decrease)
            thresh_step = thresh_step^2;
        end
        thresh = thresh * (1+thresh_step);
        last_decrease = 0; last_increase = 1;
    end
end

end