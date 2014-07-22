function nearest = MatchFDTCoordsToSurface(gm_file, inds_file, midthickness_files)
[~,~,Info,~] = BrikLoad(gm_file);
origin = Info.ORIGIN;
step_size = Info.DELTA;
inds = importdata(inds_file);
inds = inds(:,1:3);
coords = inds.*repmat(step_size,size(inds,1),1) + repmat(origin,size(inds,1),1);
coords(:,1:2) = -1*coords(:,1:2);

surf_coords = cell(2,1);
for hem = 1:2
    fid = fopen(midthickness_files{hem});
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;
    surf_coords{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(1):block_ends(1)))), 'single');
    surf_coords{hem} = permute(reshape(surf_coords{hem},[3 32492]),[2 1]);
end
surf_coords = cast(cell2mat(surf_coords),'double'); 
T = delaunayn(surf_coords);
[nearest,dists] = dsearchn(surf_coords,T,coords);

load('../data/gray.mat','orig_ind');
orig_mapping = zeros(1,2*32492);
orig_mapping(orig_ind{1}) = 1:length(orig_ind{1});
orig_mapping(orig_ind{2}+32492) = length(orig_ind{1}) + (1:length(orig_ind{2}));

nearest = orig_mapping(nearest);

end