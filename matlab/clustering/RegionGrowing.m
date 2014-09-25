function [z stats] = RegionGrowing(subject, experiment, n_clust)

stats = struct('NMI',[], 'conn_diff',[]);
loaded = load(['../data/' subject '/' experiment '.mat']);
D = loaded.D;
adj_list = loaded.adj_list;
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

instability = zeros(N,1);
for i = 1:N
    for j = adj_list{i}
        instability(i) = instability(i) + ...
                         norm(D(i, (1:N ~= i) & (1:N ~= j)) - ...
                              D(j, (1:N ~= i) & (1:N ~= j)), 2);
    end
    instability(i) = instability(i) / length(adj_list{i});
end

for i = 1:N
    instability(i) = mean([instability(i) mean(instability(adj_list{i}))]);
end

seeds = [];
parcels = {};
parcel_feat = zeros(0,N);
parcel_nbr = {};
parcel_nbr_diff = {};
added_to_parcel = zeros(N,1);
for i = 1:N
    if (all(instability(i) <= instability(adj_list{i})))
        seeds = [seeds i];
        parcels = [parcels i];
        added_to_parcel(i) = 1;
        instability(i) = -1;
        parcel_feat = [parcel_feat; D(i,:)];
        parcel_nbr = [parcel_nbr adj_list{i}];
        nbr_diff = zeros(1, length(adj_list{i}));
        for j = 1:length(adj_list{i})
            nbr_diff(j) = norm(D(i, (1:N ~= i) & (1:N ~= adj_list{i}(j))) - ...
                               D(adj_list{i}(j), (1:N ~= i) & (1:N ~= adj_list{i}(j))), 2);
        end
        parcel_nbr_diff = [parcel_nbr_diff nbr_diff];
    end
end

while(sum(added_to_parcel) < N)
    max_diff = max(cell2mat(parcel_nbr_diff));
    min_diff = min(cell2mat(parcel_nbr_diff));
    thresh = min_diff + 0.1*(max_diff - min_diff);
    
    add_mask = cell(size(parcel_nbr));
    add_vox = [];
    for i = 1:length(parcels)
        add_mask{i} = parcel_nbr_diff{i} <= thresh;
        if (sum(add_mask{i}) > 0)
            add_vox = [add_vox parcel_nbr{i}(add_mask{i})];
        end
    end
    add_vox = sort(add_vox);
    conflicts = unique(add_vox(diff(add_vox) == 0));
    
    if (~isempty(conflicts))
        for c = conflicts
            dists = cellfun(@(x, y) min([x(y == c) inf]), parcel_nbr_diff, parcel_nbr);
            [~, best_match] = min(dists);

            for i = setdiff(find(isfinite(dists)), best_match)
                add_mask{i}(parcel_nbr{i} == c) = 0;
            end
        end
    end
    
    for i = 1:length(parcels)
        % Remove neighborhood voxels added to a different parcel
        rem_vox = ismember(parcel_nbr{i}, find(added_to_parcel));
        parcel_nbr{i} = parcel_nbr{i}(~rem_vox);
        parcel_nbr_diff{i} = parcel_nbr_diff{i}(~rem_vox);
        add_mask{i} = add_mask{i}(~rem_vox);
        
        % Identify voxels to add to parcel
        add_nbr_ind = add_mask{i};
        add_nbr = parcel_nbr{i}(add_nbr_ind);
        
        % Remove voxels from neighborhood and add to parcel
        parcel_nbr{i} = parcel_nbr{i}(~add_nbr_ind);
        parcel_nbr_diff{i} = parcel_nbr_diff{i}(~add_nbr_ind);
        parcels{i} = [parcels{i} add_nbr];
        added_to_parcel(add_nbr) = 1;
        
        % Add neighbors of these voxels to neighborhood
        new_neighbors = [adj_list{add_nbr}];
        % Remove voxels already in a parcel
        new_neighbors = setdiff(new_neighbors, find(added_to_parcel));
        parcel_nbr{i} = [parcel_nbr{i} new_neighbors];
        % Remove voxels already in neighborhood
        [~,unique_ind,~] = unique(parcel_nbr{i}, 'first');
        parcel_nbr{i} = parcel_nbr{i}(sort(unique_ind));
        num_new = length(parcel_nbr{i}) - length(parcel_nbr_diff{i});
        parcel_nbr_diff{i} = [parcel_nbr_diff{i} zeros(1, num_new)];
        for j = length(parcel_nbr_diff{i}):-1:(length(parcel_nbr_diff{i}) - num_new + 1)
            parcel_nbr_diff{i}(j) = ...
                norm(parcel_feat(i,(1:N ~= seeds(i)) & (1:N ~= parcel_nbr{i}(j))) - ...
                     D(parcel_nbr{i}(j), (1:N ~= seeds(i)) & (1:N ~= parcel_nbr{i}(j))), 2);
        end
    end
end

dissim = zeros(length(parcels),length(parcels));
for m = 1:length(parcels)
    for n = 1:length(parcels)
        if (m == n)
            dissim(m,n) = inf;
        else
            dissim(m,n) = norm(parcel_feat(m, (1:N ~= seeds(m)) & (1:N ~= seeds(n))) -...
                               parcel_feat(n, (1:N ~= seeds(m)) & (1:N ~= seeds(n))), 2);
        end
    end
end

z_init = zeros(N,1);
for i = 1:length(parcels)
    z_init(parcels{i}) = i;
end

parcel_adj = inf(length(parcels),length(parcels));
for i = 1:length(parcels)
    parcel_adj(i, unique(z_init([adj_list{parcels{i}}]))) = 1;
end


parcel_vox = cellfun(@length, parcels);
while (length(parcels) > n_clust)
    [~,min_ind] = min(dissim(:).*parcel_adj(:));
    [m, n] = ind2sub(size(dissim), min_ind);
    other = (1:length(parcels) ~= m) & (1:length(parcels) ~= n);
    new_dis = (parcel_vox(m) + parcel_vox(other))./ ...
                (parcel_vox(m) + parcel_vox(n) + parcel_vox(other)) .* ...
                dissim(m,other) + ...
              (parcel_vox(n) + parcel_vox(other))./ ...
                (parcel_vox(m) + parcel_vox(n) + parcel_vox(other)) .* ...
                dissim(n,other) + ...
              - parcel_vox(other)./ ...
                (parcel_vox(m) + parcel_vox(n) + parcel_vox(other)) .* ...
                dissim(m,n);
    dissim(m,other) = new_dis;
    dissim(other,m) = new_dis;
    dissim = dissim(1:length(parcels) ~= n,1:length(parcels) ~= n);
    
    parcel_adj(m,parcel_adj(n,:)==1) = 1;
    parcel_adj = parcel_adj(1:length(parcels) ~= n,1:length(parcels) ~= n);
    
    parcels{m} = [parcels{m} parcels{n}];
    parcels = parcels(1:length(parcels) ~= n);
end

z = zeros(N,1);
for i = 1:n_clust
    z(parcels{i}) = i;
end

if (~isempty(gt_z))
    stats.NMI = CalcNMI(gt_z, z);
end
end