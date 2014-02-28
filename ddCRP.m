function [map_z stats] = ddCRP(D, adj_list, coords, ...
                               init_c, const_c, labels, bold, gt_z, ...
                               num_passes, alpha, kappa, nu, sigsq, ...
                               stats_interval, verbose)

hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq);
pid = randi(10000);
nvox = size(coords, 1);

if (isempty(const_c))
    const_c = zeros(nvox, 1);
end

if (isempty(init_c))
    init_c = zeros(nvox, 1);
end

c = const_c;
for i = find(const_c==0)'
    if (init_c(i) == 0)
        c(i) = adj_list{i}(randi(length(adj_list{i})));
    else
        c(i) = init_c(i);
    end
end

G = sparse(1:nvox,c,1,nvox,nvox);
[K, z, parcels] = ConnectedComp(G);
curr_lp = FullProbabilityddCRP(D, c, parcels, alpha, hyp);


stats = struct('times',[],'lp',[],'NMI',[],'K',[], ...
               'conn_diff', zeros(0,4), 'z', zeros(0,nvox));
max_lp = -Inf;
t0 = cputime;
steps = 0;
for pass = 1:num_passes
    nonconst_vox = find(const_c==0);
    order = nonconst_vox(randperm(length(nonconst_vox)))';
    
    for i = order
        if (curr_lp > max_lp)
            max_lp = curr_lp;
            map_z = z;
        end
        
        if (mod(steps, stats_interval) == 0)
            stats = UpdateStats(stats, t0, curr_lp, K, z, steps, gt_z, map_z, pid, verbose);
        end
        
        if (~isempty(labels) && ~isempty(bold))
            stats.conn_diff = [stats.conn_diff; ...
                               CalcPPAConnDiff(z, labels, coords, bold)];
        end
        
        if (c(i) == i)
            rem_delta_lp = -log(alpha);
            z_rem = z; parcels_rem = parcels;
        else
            G(i,c(i)) = 0;
            [K_rem, z_rem, parcels_rem] = ConnectedComp(G);
            if (K_rem ~= K)
                rem_delta_lp = -LikelihoodDiff(D, ...
                                  parcels_rem, z_rem(i), z_rem(c(i)), hyp);
            else
                rem_delta_lp = 0;
            end
        end
        
        adj_list_i = adj_list{i};
        lp = zeros(length(adj_list_i)+1, 1);
        lp(end) = log(alpha);
        for n_ind = 1:length(adj_list_i)
            n = adj_list_i(n_ind);
            if (z_rem(n) == z_rem(c(i)))  % Clustered with old neighbor
                lp(n_ind) = -rem_delta_lp;
            elseif (z_rem(n) ~= z_rem(i))  % Not already clustered with n
                lp(n_ind) = LikelihoodDiff(D, parcels_rem, z_rem(i), z_rem(n), hyp);
            end
        end
        
        new_neighbor = ChooseFromLP(lp);
        if (new_neighbor <= length(adj_list_i))
            c(i) = adj_list_i(new_neighbor);
        else
            c(i) = i;
        end
        curr_lp = curr_lp + rem_delta_lp + lp(new_neighbor);
        G(i,c(i)) = 1;
        [K, z, parcels] = ConnectedComp(G);
        steps = steps + 1;
    end
end

stats = UpdateStats(stats, t0, curr_lp, K, z, steps, gt_z, map_z, pid, verbose);

end

function [K, z, parcels] = ConnectedComp(G)
    [K, z] = graphconncomp(G, 'Weak', true);
    [sorted_z, sorted_i] = sort(z);
    parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (K+1)]))));
end

function ld = LikelihoodDiff(D, parcels_split, split_i1, split_i2, hyp)
    K = length(parcels_split);
    s = zeros(2*K, 3);
    for i = 1:K
        samples = D(parcels_split{i}, parcels_split{split_i1});
        if (i == split_i1)
            samples = samples(logical(triu(ones(size(samples)),1)));
            if (isempty(samples))
                continue;
            end
        else
            samples = samples(:);
        end
        s(i,1) = length(samples);
        s(i,2) = sum(samples)/s(i,1);
        s(i,3) = sum((samples-s(i,2)).^2);
    end
    for i = 1:K
        samples = D(parcels_split{i}, parcels_split{split_i2});
        if (i == split_i1)
            continue;
        end
        if (i == split_i2)
            samples = samples(logical(triu(ones(size(samples)),1)));
            if (isempty(samples))
                continue;
            end
        else
            samples = samples(:);
        end
        s(i+K,1) = length(samples);
        s(i+K,2) = sum(samples)/s(i+K,1);
        s(i+K,3) = sum((samples-s(i+K,2)).^2);
    end
    split_ll = LogLikelihood(s, hyp);
    
    m = zeros(K, 3);
    for i = 1:K
        if (i == split_i1)
            i11 = i;
            i12 = split_i2;
            i22 = split_i2 + K;
            m(i,1) = s(i11,1) + s(i12,1) + s(i22,1);
            m(i,2) = (s(i11,1)*s(i11,2) + s(i12,1)*s(i12,2) + s(i22,1)*s(i22,2))/m(i,1);
            if (s(i11,1) + s(i22,1) > 0)
                %Combine i11 and i22
                m(i,3) = s(i11,3) + s(i22,3) + (s(i11,1)*s(i22,1))/(s(i11,1)+s(i22,1)) * (s(i11,2)- s(i22,2))^2;
                mean_11_22 = (s(i11,1)*s(i11,2)+s(i22,1)*s(i22,2))/(s(i11,1)+s(i22,1));
                %Combine with i12
                m(i,3) = m(i,3) + s(i12,3) + (s(i12,1)*(s(i11,1)+s(i22,1)))/m(i,1) * (s(i12,2)- mean_11_22)^2;
            else
                m(i,3) = 0;
            end
        elseif (i ~= split_i2)
            m(i,1) = s(i,1) + s(i+K,1);
            m(i,2) = (s(i,1)*s(i,2) + s(i+K,1)*s(i+K,2))/m(i,1);
            m(i,3) = s(i,3) + s(i+K,3) + (s(i,1)*s(i+K,1))/m(i,1) * (s(i,2)- s(i+K,2))^2;
        end
    end
    merge_ll = LogLikelihood(m, hyp);
    
    ld = merge_ll - split_ll;
end

function stats = UpdateStats(stats, t0, curr_lp, K, z, steps, gt_z, map_z, pid, verbose)
    stats.times = [stats.times (cputime-t0)];
    stats.lp = [stats.lp curr_lp];
    stats.K = [stats.K K];
    stats.z = [stats.z; z];
    if (verbose)
        disp(['Step: ' num2str(steps) ...
              '  Time: ' num2str(cputime-t0) ...
              '  LP: ' num2str(curr_lp)]);
    end
    if (~isempty(gt_z))
        stats.NMI = [stats.NMI CalcNMI(gt_z, map_z)];
    end
    save(['/data/supervoxel/output/temp/' num2str(pid) '.mat'], ...
        'map_z', 'stats');
end