function [z c] = GenerateRandomClustering( adj_list, K, z_init)

    N = length(adj_list);
    if (nargin < 3)
        merge_vox = 1:N;
        K_curr = N;
        z = 1:N;
    else
        z = z_init;
        merge_vox = find(z_init==0);
        K_curr = length(merge_vox);
        z(merge_vox) = (max(z_init)+1):(max(z_init)+K_curr);
    end
    while (K_curr > K)
%         if (mod(K_curr, 100) == 0)
%             disp(K_curr);
%         end
        clustLabels = unique(z(merge_vox));
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

