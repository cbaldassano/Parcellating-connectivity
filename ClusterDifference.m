function diff_map = ClusterDifference(z1, z2)

z1_pairs = PairMat(z1);
z2_pairs = PairMat(z2);

diff_map = sum(xor(z1_pairs, z2_pairs));

end

function z_pairs = PairMat(z)
z_pairs = false(length(z));

[sorted_z, sorted_i] = sort(z);
bins = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

for c = 1:length(bins)
    repbin = repmat(bins{c}, length(bins{c}), 1);
    z_pairs(sub2ind([length(z) length(z)], ...
        reshape(repbin, 1, size(repbin,1)^2), ...
        reshape(repbin', 1, size(repbin,1)^2))) = true;
end

end