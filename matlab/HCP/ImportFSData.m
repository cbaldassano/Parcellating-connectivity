function map = ImportFSData(pial_files, surf_data_files, coords)
%ImportFSData({'../data/wmRetinotopy/fsaverage.L.pial.surf.gii', '../data/wmRetinotopy/fsaverage.R.pial.surf.gii'}, {'../data/wmRetinotopy/ecc-lh-fsavg-sm2.gii', '../data/wmRetinotopy/ecc-rh-fsavg-sm2.gii'}, coords );
data_L = LoadMetric(surf_data_files{1});
data_R = LoadMetric(surf_data_files{2});

fscoords_L = reshape(LoadMetric(pial_files{1}),3,163842)';
fscoords_R = reshape(LoadMetric(pial_files{2}),3,163842)';

T = delaunayn(cast(coords,'double'));
nearest_L = dsearchn(cast(coords,'double'),T,cast(fscoords_L,'double'));
nearest_R = dsearchn(cast(coords,'double'),T,cast(fscoords_R,'double'));

map = zeros(59412,1);
for v = 1:length(map)
    if (mod(v,1000) == 1)
        fprintf('%d ', v);
    end
    mapped_data_L = find(nearest_L == v);
    mapped_data_R = find(nearest_R == v);
    if (~isempty(mapped_data_L) || ~isempty(mapped_data_R))
        map(v) = mean([data_L(mapped_data_L); data_R(mapped_data_R)]);
    end
end
end