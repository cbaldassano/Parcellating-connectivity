function degree_density = MeanDegree(Dfile, map_z, parcelnums)

full = matfile(Dfile);

degree_density = zeros(length(parcelnums),1);
for i = 1:length(parcelnums)
    disp(['Calculating density for ' num2str(parcelnums(i))]);
    parcel_inds = find(map_z==parcelnums(i));
    nonparcel_inds = map_z~=parcelnums(i);
    for j = parcel_inds
        vox_conn = full.D(j,:);
        degree_density(i) = degree_density(i) + ...
            sum(vox_conn(nonparcel_inds))/(sum(nonparcel_inds)*length(parcel_inds));
    end
end

end