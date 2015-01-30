function [total_conn_LR num_vox_LR] = DegreeDist(Dfile)
% For anat: '/data/supervoxel/data/Q3/group/full.mat'

load('/data/supervoxel/data/group468/coords.mat');
load('../../output/group468/ddCRP3000_rng1.mat');
parcelnums = [55 54 60 61 13 8 17 16; ...
               141 140 145 146 98 93 101 100]';

Nbins = 25;
full = matfile(Dfile);
dist_edges = linspace(0,250,Nbins+1);

total_conn = zeros(length(parcelnums(:)),Nbins);
num_vox = zeros(length(parcelnums(:)),Nbins);
for i = 1:length(parcelnums(:))
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
        total_conn(i,:) = total_conn(i,:) + cellfun(@(x) sum(vox_conn(x)), bin_inds);
        num_vox(i,:) = num_vox(i,:) + cellfun(@length, bin_inds);
    end
end

total_conn_LR = zeros(8, Nbins);
num_vox_LR = zeros(8, Nbins);
for i = 1:8
   total_conn_LR(i,:) =  (total_conn(i,:) + total_conn(i+8,:))/sum(map_z == parcelnums(i,1) | map_z == parcelnums(i,2));
   num_vox_LR(i,:) = (num_vox(i,:) + num_vox(i+8,:))/sum(map_z == parcelnums(i,1) | map_z == parcelnums(i,2));
end

end
