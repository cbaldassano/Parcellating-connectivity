function ShowPairProbs(map_z, pair_probs)

pp = squareform(pair_probs);
K = max(map_z);
N = length(map_z);

clust_col = repmat(PTPalette(12),ceil(K/12),1);
clust_col = clust_col(1:K,:);

[sorted_z, sorted_i] = sort(map_z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (K+1)]))));

vox_col = zeros(N,3);
for i = 1:N
    clust_mean = cellfun(@(x) mean(pp(i,x)), parcels);
    clust_mean = clust_mean/sum(clust_mean);
    vox_col(i,:) = clust_mean*clust_col;
end

image(reshape(vox_col,18,18,3));
end