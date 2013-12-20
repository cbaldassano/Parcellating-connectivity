function DivideSurface(midthickness_files, surface_template_file, z, orig_ind, overlay, save_prefix, overlay_suffix)

[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));
orig_parcels = cell(2,1);
for i = 1:length(parcels)
    if (parcels{i}(1) <= 29696)
        orig_parcels{1} = [orig_parcels{1} {orig_ind{1}(parcels{i})}];
    else
        orig_parcels{2} = [orig_parcels{2} {orig_ind{2}(parcels{i} - 29696)}];
    end
end

max_overlay = max(max(overlay{1}),max(overlay{2}));
colors = rand(max_overlay + 2, 3);
colors(1,:) = [0 0 0];
colors(2,:) = [1 1 1];

% Inspired by Guillaume Flandin <Guillaume@artefact.tk> http://www.artefact.tk/software/matlab/gifti/
fid = fopen(surface_template_file);
template_xml = fread(fid,'*char')';
fclose(fid);

fid_filelist = fopen([save_prefix '_flist.txt'], 'w');
fid_sizelist = fopen([save_prefix '_slist.txt'], 'w');
coords = cell(2,1);
tri = cell(2,1);
for hem = 1:2
    disp(['Loading surface ' num2str(hem) '...']);
    fid = fopen(midthickness_files{hem});
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    
    disp('   Extracting data...');
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;
    
    coords{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(1):block_ends(1)))), 'single');
    
    coords{hem} = permute(reshape(coords{hem},[3 32492]),[2 1]);
    
    tri{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(2):block_ends(2)))), 'int32');
    
    tri{hem} = permute(reshape(tri{hem},[3 64980]),[2 1]) + 1;
    
    for i = 1:length(orig_parcels{hem})
        parcel_tri = tri{hem}(any(ismember(tri{hem}, orig_parcels{hem}{i}), 2), :);
        if (isempty(parcel_tri))
            continue;
        end
        parcel_overlay = overlay{hem}(parcel_tri);
        parcel_overlay(~ismember(parcel_tri, orig_parcels{hem}{i})) = -1;
        fid = fopen([save_prefix '_' num2str(hem) '_' num2str(i) '_' overlay_suffix '.clrs'], 'w');
        for f = 1:size(parcel_overlay,1)
            for v = 1:3
                val = parcel_overlay(f,v) + 2;
                fprintf(fid, '%f %f %f\n', colors(val,1), colors(val,2), colors(val,3));
            end
        end
        fclose(fid);
        fprintf(fid_sizelist, '%d, ', 9*size(parcel_overlay,1));
        [~,filename] = fileparts([save_prefix '_' num2str(hem) '_' num2str(i)]);
        fprintf(fid_filelist, '''%s'', ', filename);
        
        verts = unique(parcel_tri(:));
        inverse_ind = zeros(32492,1);
        inverse_ind(verts) = 1:length(verts);
        parcel_tri = cast(permute(inverse_ind(parcel_tri) - 1, [2 1]), 'int32');
        parcel_coords = permute(coords{hem}(verts, :), [2 1]);
        
        fid = fopen([save_prefix '_' num2str(hem) '_' num2str(i) '.gii'], 'w');
        fprintf(fid, template_xml, ...
            size(parcel_coords,2), base64encode(dzip(typecast(parcel_coords(:), 'uint8'))), ...
            size(parcel_tri,2), base64encode(dzip(typecast(parcel_tri(:), 'uint8'))));
        fclose(fid);
        system(['mris_convert ' save_prefix '_' num2str(hem) '_' num2str(i) '.gii ' save_prefix '_' num2str(hem) '_' num2str(i) '.fsm']);
    end
end

fclose(fid_filelist);
fclose(fid_sizelist);
end

