function WriteTRK(pts, mri_info_file, savename)
% pts is a cell array of m x 3 arrays giving ordered points for each track
fid_info = fopen(mri_info_file);
dims = textscan(fid_info, '%s %d x %d x %d x %d', 1, 'HeaderLines', 2);
dims = cell2mat(dims([2 3 4]));
voxel_sizes = textscan(fid_info, '%s %s %f, %f, %f', 1);
voxel_sizes = cell2mat(voxel_sizes(3:end));
orient = textscan(fid_info, '%s %s', 1, 'HeaderLines', 14, 'Delimiter', ':');
orient = orient{2}{1};
vox_to_ras = textscan(fid_info, '%f %f %f %f', 4, 'HeaderLines', 4);
vox_to_ras = cell2mat(vox_to_ras);
fclose(fid_info);

fid = fopen(savename, 'w');

% id_string
fwrite(fid, 'TRACKM', 'char*1');

%dim
fwrite(fid, dims, 'int16');

% voxel_size
fwrite(fid, voxel_sizes, 'float');

% origin
fwrite(fid, [0 0 0], 'float');

% n_scalars
fwrite(fid, 0, 'int16');

% scalar_name
fwrite(fid, zeros(10,20), 'char*1');

% n_properties
fwrite(fid, 0, 'int16');

% property_name
fwrite(fid, zeros(10,20), 'char*1');

% vox_to_ras
fwrite(fid, vox_to_ras', 'float');

% reserved
fwrite(fid, zeros(444,1), 'char*1');

% voxel_order
fwrite(fid, orient, 'char*1'); fwrite(fid, zeros(1,1), 'char*1');

% pad2
fwrite(fid, zeros(4,1), 'char*1');

% image_orientation_patient
fwrite(fid, zeros(6,1), 'float');

% pad1
fwrite(fid, zeros(2,1), 'char*1');

% invert/swap
fwrite(fid, zeros(6,1), 'char*1');

% n_count
fwrite(fid, length(pts), 'int');

% version
fwrite(fid, 2, 'int');

%hdr_size
fwrite(fid, 1000, 'int');

for i = 1:length(pts)
    fwrite(fid, size(pts{i},1), 'int');
    for j = 1:size(pts{i},1)
        fwrite(fid, pts{i}(j,:), 'float');
    end
end
fclose(fid);
end