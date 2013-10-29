function DivideSurface(midthickness_files, surface_template_file, z, orig_ind, save_prefix)

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

% Inspired by Guillaume Flandin <Guillaume@artefact.tk> http://www.artefact.tk/software/matlab/gifti/
fid = fopen(surface_template_file);
template_xml = fread(fid,'*char')';
fclose(fid);

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

system(['ls -1v ' save_prefix '*.fsm | sed s/^.\*\\/\// > ' save_prefix '.txt']);

end

