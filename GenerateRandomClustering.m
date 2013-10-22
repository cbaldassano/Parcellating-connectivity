function [z c] = GenerateRandomClustering( adj_list, K )

    N = length(adj_list);
    z = 1:N;
    K_curr = N;
    while (K_curr > K)
        if (mod(K_curr, 1000) == 0)
            disp(K_curr);
        end
        clustLabels = unique(z);
        zToMerge = clustLabels(randi(length(clustLabels),1));
        mergeChoices = unique(z([adj_list{z==zToMerge}]));
        mergeChoices = setdiff(mergeChoices,zToMerge);
        if (~isempty(mergeChoices))
            z(z==zToMerge) = mergeChoices(randi(length(mergeChoices),1));
            K_curr = K_curr - 1;
        end
    end

    zRelabel = zeros(N,1);
    relabelInd = 1;
    for v=unique(z)
        zRelabel(z==v) = relabelInd;
        relabelInd = relabelInd+1;
    end
    z = zRelabel;
    if (nargout == 2)
        c = ClusterSpanningTrees(z, adj_list);
    end
end

