function mean_sq_diff = CalcMeanNeighborDiff(D, adj_list)

neighbor_diffs = zeros(length(adj_list),1);
for i = 1:length(adj_list)
    diff_vecs = (cast(D(adj_list{i},:),'double') - repmat(cast(D(i,:),'double'),length(adj_list{i}),1)).^2;
    diff_vecs(isinf(diff_vecs)) = 0;
    ref_vec = cast(D(i,:),'double').^2;
    ref_vec(isinf(ref_vec)) = 0;
    neighbor_diffs(i) = mean(sum(diff_vecs, 2))/sum(ref_vec);
end

mean_sq_diff = mean(neighbor_diffs);
end