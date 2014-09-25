function c = ClusterSpanningTrees(z, adj_list)

    nvox = length(adj_list);
    for i = 1:nvox
        adj_list{i} = adj_list{i}(z(adj_list{i}) == z(i));
        adj_list{i} = adj_list{i}(randperm(length(adj_list{i})));
    end
    neighbor_count = cellfun(@length, adj_list);
    node_list = zeros(sum(neighbor_count), 1);
    next_edge = 1;
    for i = 1:nvox
        if (neighbor_count(i) > 0)
            node_list(next_edge:(next_edge+neighbor_count(i)-1)) = i;
            next_edge = next_edge + neighbor_count(i);
        end
    end
    G = sparse(node_list, [adj_list{:}]', 1, nvox, nvox);
    
    c = zeros(length(adj_list),1);
    for clust = unique(z)';
        clust_vox = find(z==clust);
        [~,parents] = graphminspantree(G,clust_vox(randi(length(clust_vox),1)));
        c(clust_vox) = parents(clust_vox);
    end
    roots = find(c==0);
    c(roots) = roots;

end

