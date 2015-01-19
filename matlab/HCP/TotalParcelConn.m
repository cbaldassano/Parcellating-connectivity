function conn = TotalParcelConn(Dfile, map_z, parcelnums)

full = matfile(Dfile);

[sorted_bins, sorted_i] = sort(map_z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_bins (max(map_z)+1)]))));
conn = zeros(length(parcelnums), length(parcels));

for i = 1:length(parcelnums)
    parcel_inds = find(map_z==parcelnums(i));
    for j = parcel_inds
        vox_conn = full.D(j,:);
        conn(i,:) = conn(i,:) + cellfun(@(x) mean(vox_conn(x)), parcels)/length(parcel_inds);
    end
end

end