function overlays = AssignToSurface(midthickness_files, maps)

[~,~,Info] = BrikLoad(maps{1});
V_dict = cell(2,1);
origin = Info.ORIGIN;
step_size = Info.DELTA;
coords = cell(2,1);
for hem = 1:2
    fid = fopen(midthickness_files{hem});
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;
    coords{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(1):block_ends(1)))), 'single');
    coords{hem} = permute(reshape(coords{hem},[3 32492]),[2 1]);
    coords{hem}(:,1:2) = -1*coords{hem}(:,1:2);
    
    V_dict{hem} = zeros(Info.DATASET_DIMENSIONS(1:3));
    for v = 1:size(coords{hem},1)
        V_ind = round((coords{hem}(v,:)-origin)./step_size) + 1;
        V_dict{hem}(V_ind(1), V_ind(2), V_ind(3)) = v;
    end
end

overlays = cell(length(maps),1);
for m = 1:length(maps)
    [~,V,~,~] = BrikLoad(maps{m});
    ind = find(V > 1);
    
    overlays{m} = cell(2,1);
    for hem = 1:2
        overlays{m}{hem} = zeros(32492, 1);
        surf_ind = V_dict{hem}(ind);
        surf_ind = surf_ind(surf_ind > 0);
        overlays{m}{hem}(surf_ind) = 1;
    end
end

end

