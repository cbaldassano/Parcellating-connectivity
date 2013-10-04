function [map_z stats] = ddCRP(subject, experiment, num_passes, ...
                               alpha, kappa, nu, sigsq)

hyp = [kappa nu sigsq];
stats = struct('times',[],'lp',[],'NMI',[],'K',[], 'conn_diff',zeros(0,4));

loaded = load(['../data/' subject '/' experiment '.mat']);
D = loaded.D;
adj_list = loaded.adj_list;
coords = loaded.coords;
nvox = size(coords, 1);

if (strcmp(experiment, 'PPA'))
    const_c = loaded.const_c;
    labels = loaded.labels;
    bold = loaded.bold;
else
    const_c = zeros(nvox, 1);
end
if (strcmp(subject, 'synth'))
    gt_z = loaded.z;
else
    gt_z = [];
end
clear loaded;


c = const_c;
for i = find(const_c==0)'
    c(i) = adj_list{i}(randi(length(adj_list{i})));
end

G = sparse(1:nvox,c,1,nvox,nvox);
[K, z, parcels] = ConnectedComp(G);
curr_lp = FullProbabilityddCRP(D, c, parcels, alpha, hyp);

max_lp = -Inf;
tic;
for pass = 1:num_passes
    nonconst_vox = find(const_c==0);
    order = nonconst_vox(randperm(length(nonconst_vox)))';
    
    for i = order
        if (curr_lp > max_lp)
            max_lp = curr_lp;
            map_z = z;
        end
        
        stats.times = [stats.times toc];
        stats.lp = [stats.lp curr_lp];
        stats.K = [stats.K K];
        if (~isempty(gt_z))
            stats.NMI = [stats.NMI CalcNMI(gt_z,MAPz)];
            if (abs(stats.NMI(end)-1)<10^(-8))
                return;
            end
        end
        
        if (strcmp(experiment,'PPA'))
            stats.conn_diff = [conn_diff; ...
                               CalcPPAConnDiff(z, labels, coords, bold)];
        end
        
        if (c(i) == i)
            rem_delta_lp = -log(alpha);
        else
            rem_delta_lp = 0;
        end
        G(i,c(i)) = 0;
        [K_rem, z_rem, parcels_rem] = ConnectedComp(G);
        if (K_rem ~= K)
            rem_delta_lp = rem_delta_lp - ...
                      LikelihoodDiff(parcels_rem, z_rem(i), z_rem(c(i)),...
                                     parcels, z(i), hyp);
        end
        
        lp = zeros(length(adj_list{i})+1, 1);
        lp(end) = log(alpha);
        for n_ind = 1:length(adj_list{i})
            n = adj_list{i}(n_ind);
            c(i) = n;
            G_n = sparse(1:N,c,1,N,N);
            [K_n, z_n, parcels_n] = ConnectedComp(G_n);
            if (K_n == K_rem)
                continue;
            end
            
            lp(n_ind) = LikelihoodDiff(parcels_rem, z_rem(i), z_rem(n), ...
                                       parcels_n, z_n(i), hyp);
        end
        
        new_neighbor = ChooseFromLP(lp);
        if (new_neighbor <= length(adj_list{i}))
            c(i) = adj_list{i}(new_neighbor);
        else
            c(i) = i;
        end
        curr_lp = curr_lp + rem_delta_lp + lp(new_neighbor);
        G(i,c(i)) = 1;
        [K, z, parcels] = ConnectedComp(G);
    end
end
end

function [K, z, parcels] = ConnectedComp(G)
    [K, z] = graphconncomp(G, 'Weak', true);
    parcels = arrayfun(@(x) find(z==x), 1:K, 'UniformOutput', false);
end

function ld = LikelihoodDiff(parcels_split, split_i1, split_i2, ...
                             parcels_merge, merge_i, hyp)
    split_ll = 0;
    for i = 1:length(parcels_split)
        split_ll = split_ll + ...
            LogLikelihoodOfClusterPair(D, parcels_split{i}, ...
                                       parcels_split{split_i1}, hyp);
        if (i ~= split_i1)
            split_ll = split_ll + ...
                LogLikelihoodOfClusterPair(D, parcels_split{i}, ...
                                           parcels_split{split_i2}, hyp);
        end
    end
    
    merge_ll = 0;
    for i = 1:length(parcels_merge)
        merge_ll = merge_ll + ...
            LogLikelihoodOfClusterPair(D, parcels_merge{i}, ...
                                       parcels_merge{merge_i}, hyp);
    end
    
    ld = merge_ll - split_ll;
end

function ll = LogLikelihoodOfClusterPair(D, inds1, inds2, hyp)
    samples = D(inds1,inds2);
    if (length(inds1)==length(inds2) && all(inds1==inds2))
            samples = samples(logical(triu(ones(size(samples)),1)));
    else
            samples = samples(:);
    end
    
    ll = LogLikelihood(samples, hyp);
end
