function [z c] = GenerateRandomClustering( adj_list, K )

    N = length(adj_list);
    z = 1:N;
    while (length(unique(z)) > K)
        clustLabels = unique(z);
        zToMerge = clustLabels(randi(length(clustLabels),1));
        mergeChoices = unique(z([adj_list{z==zToMerge}]));
        mergeChoices = setdiff(mergeChoices,zToMerge);
        if (~isempty(mergeChoices))
            z(z==zToMerge) = mergeChoices(randi(length(mergeChoices),1));
        end
    end

    zRelabel = zeros(N,1);
    relabelInd = 1;
    for v=unique(z)
        zRelabel(z==v) = relabelInd;
        relabelInd = relabelInd+1;
    end
    z = zRelabel;
    c = ClusterSpanningTrees(z, adj_list);
end

