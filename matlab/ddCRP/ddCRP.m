function [map_z stats pair_prob] = ddCRP(D, adj_list, ...
                                    init_c, const_c, gt_z, ...
                                    num_passes, alpha, kappa, nu, sigsq, ...
                                    burn_in_passes, ...
                                    stats_interval, verbose)

hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq);
pid = randi(10000);
nvox = length(adj_list);

if (isempty(const_c))
    const_c = zeros(nvox, 1);
end

if (isempty(init_c))
    init_c = zeros(nvox, 1);
end

if (~isempty(burn_in_passes))
    pair_prob = zeros(1,nvox*(nvox-1)/2);
end

c = const_c;
for i = find(const_c==0)'
    if (init_c(i) == 0)
        neighbors = [adj_list{i} i];
        c(i) = neighbors(randi(length(neighbors)));
    else
        c(i) = init_c(i);
    end
end

G = sparse(1:nvox,c,1,nvox,nvox);
[K, z, parcels] = ConnectedComp(G);

sym = CheckSymApprox(D);
      
curr_lp = FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym);


stats = struct('times',[],'lp',[],'NMI',[],'K',[], 'z', zeros(0,nvox), 'c', zeros(0,nvox));
max_lp = -Inf;
t0 = tic;
steps = 0;
for pass = 1:num_passes
    nonconst_vox = find(const_c==0);
    order = nonconst_vox(randperm(length(nonconst_vox)))';
    
    for i = order
        if (curr_lp > max_lp)
            max_lp = curr_lp;
            map_z = z;
        end
        
        if (~isempty(burn_in_passes) && pass > burn_in_passes)
            pair_prob = pair_prob + (1 - pdist(z', 'hamming'));
        end
        
        if (mod(steps, stats_interval) == 0)
            stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, pid, verbose);
        end
        
        G(i,c(i)) = 0;
        if (c(i) == i)
            rem_delta_lp = -log(alpha);
            z_rem = z; parcels_rem = parcels;
        else
            [K_rem, z_rem, parcels_rem] = ConnectedComp(G);
            if (K_rem ~= K)
                rem_delta_lp = -LikelihoodDiff(D, ...
                                  parcels_rem, z_rem(i), z_rem(c(i)), hyp, sym);
            else
                rem_delta_lp = 0;
            end
        end
        
        adj_list_i = adj_list{i};
        lp = zeros(length(adj_list_i)+1, 1);
        lp(end) = log(alpha);
        cached_merge = zeros(length(adj_list_i),1);
        for n_ind = 1:length(adj_list_i)
            n = adj_list_i(n_ind);
            if (z_rem(n) == z_rem(c(i)))
                % Just undoing edge removal
                lp(n_ind) = -rem_delta_lp - (c(i)==i)*log(alpha);
            elseif (z_rem(n) ~= z_rem(i)) 
                % Proposing novel merge
                % First check cache to see if this is already computed
                prev_lp = find(cached_merge == z_rem(n),1);
                if (~isempty(prev_lp))
                    lp(n_ind) = lp(prev_lp);
                else
                    lp(n_ind) = LikelihoodDiff(D, parcels_rem, z_rem(i), z_rem(n), hyp, sym);
                    cached_merge(n_ind) = z_rem(n);
                end
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

stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, pid, verbose);
if (~isempty(burn_in_passes))
    pair_prob = pair_prob / (sum(const_c==0)*(num_passes-burn_in_passes));
else
    pair_prob = [];
end

end

function [K, z, parcels] = ConnectedComp(G)
    [K, z] = graphconncomp(G, 'Weak', true);
    [sorted_z, sorted_i] = sort(z);
    parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (K+1)]))));
end

function ld = LikelihoodDiff(D, parcels_split, split_i1, split_i2, hyp, sym)
    K = length(parcels_split);
    s = zeros(K, K, 3);
    for split_ind = [split_i1 split_i2]
        for i = 1:K
            samples = D(parcels_split{i}, parcels_split{split_ind});
            if (sym && i == split_ind)
                samples = samples(logical(triu(ones(size(samples)),1)));
            else
                samples = samples(:);
            end
            s(i,split_ind,:) = SufficientStats(samples);
        end
        if (~sym)
            for i = 1:K
                samples = D(parcels_split{split_ind}, parcels_split{i});
                if (i == split_ind)
                    off_diags = true(size(samples));
                    off_diags(1:(size(samples,1)+1):end) = false;
                    samples = samples(off_diags);
                else
                    samples = samples(:);
                end
                s(split_ind,i,:) = SufficientStats(samples);
            end
        end
    end
    
    if (sym)
        split_ll = LogLikelihood([...
            reshape(s(:,split_i1,:),[],3); ...
            reshape(s(1:K ~= split_i1, split_i2,:),[],3)], hyp);
    else
        split_ll = LogLikelihood([...
            reshape(s(:, split_i1,:),[],3); ...
            reshape(s(:, split_i2,:),[],3); ...
            reshape(s(split_i1, (1:K ~= split_i1) & (1:K ~= split_i2),:),[],3); ...
            reshape(s(split_i2, (1:K ~= split_i1) & (1:K ~= split_i2),:),[],3)], hyp);
    end
    
    m = zeros(2, K, 3);
    for dir = 1:2
        if (dir == 2)
            if (sym)
                break;
            else
                s = permute(s, [2 1 3]);
            end
        end
        
        for i = 1:K
            if (i ~= split_i1 && i ~= split_i2)
                s_m = reshape(s(i, [split_i1 split_i2], :),2,3);
                m(dir,i,:) = MergeSuffStats(s_m);
                
            end
        end
    end
    if (sym)
        m_central = MergeSuffStats(...
            [MergeSuffStats(reshape(s(split_i1, [split_i1 split_i2], :),2,3)); ...
            reshape(s(split_i2, split_i2, :),1,3)]);
        merge_ll = LogLikelihood([reshape(m(1, :, :),[],3); m_central], hyp);
    else
        m_central = MergeSuffStats(...
            [MergeSuffStats(reshape(s(split_i1, [split_i1 split_i2], :),2,3)); ...
             MergeSuffStats(reshape(s(split_i2, [split_i1 split_i2], :),2,3))]);
         merge_ll = LogLikelihood([reshape(m(1, :, :),[],3); reshape(m(2, :, :),[],3); m_central], hyp);
    end
    
    ld = merge_ll - split_ll;
end

function suffstats = SufficientStats(samples)
    suffstats = zeros(3,1);
    if (isempty(samples))
        return;
    end
    suffstats(1) = length(samples);
    suffstats(2) = sum(samples)/suffstats(1);
    suffstats(3) = sum((samples-suffstats(2)).^2);
end

function m = MergeSuffStats(s_m)
    m = zeros(1,3);
    m(1) = s_m(1,1) + s_m(2,1);
    m(2) = (s_m(1,1)*s_m(1,2) + s_m(2,1)*s_m(2,2))/m(1);
    m(3) = s_m(1,3) + s_m(2,3) + ...
                 (s_m(1,1)*s_m(2,1))/m(1) * (s_m(1,2) - s_m(2,2))^2;
end

function stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, pid, verbose)
    stats.lp = [stats.lp curr_lp];
    stats.K = [stats.K K];
    stats.z = [stats.z; z];
    elapsed = toc(t0);
    stats.times = [stats.times elapsed];
    stats.c = [stats.c; c'];
    if (verbose)
        disp(['Step: ' num2str(steps) ...
              '  Time: ' num2str(elapsed) ...
              '  LP: ' num2str(curr_lp) ...
              '  K: ' num2str(K)]);
    end
    if (~isempty(gt_z))
        stats.NMI = [stats.NMI CalcNMI(gt_z, map_z)];
    end
    save(['/data/supervoxel/output/temp/' num2str(pid) '.mat'], ...
        'map_z', 'stats');
end