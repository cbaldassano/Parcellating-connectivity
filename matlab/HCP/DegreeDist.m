function degree_density = DegreeDist(Dfile, coords, map_z, parcelnums)

Nbins = 25;
full = matfile(Dfile);
dist_edges = linspace(0,250,Nbins+1);

degree_density = zeros(length(parcelnums),Nbins);
for i = 1:length(parcelnums)
    disp(['Calculating density for ' num2str(parcelnums(i))]);
    parcel_inds = find(map_z==parcelnums(i));
    nonparcel_inds = map_z~=parcelnums(i);
    for j = parcel_inds
        vox_conn = full.D(j,:);
        vox_conn = vox_conn(nonparcel_inds);
        coord_dists = sqrt(sum((bsxfun(@minus, coords(nonparcel_inds,:), coords(j,:))).^2,2));
        [~, coord_bins] = histc(coord_dists, dist_edges);
        [sorted_bins, sorted_i] = sort(coord_bins');
        bin_inds = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_bins (Nbins+1)]))));
        bin_inds = [bin_inds cell(1,Nbins-length(bin_inds))];
        degree_density(i,:) = degree_density(i,:) + ...
            cellfun(@(x) sum(vox_conn(x)), bin_inds)/(sum(nonparcel_inds)*length(parcel_inds));
    end
end

end