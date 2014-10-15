function overlays = ClosestSurfaceVertex(midthickness_files, maps, max_dist, thresh)

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
    
end
T = delaunayn(double(cell2mat(coords)));

overlays = cell(length(maps),1);
for m = 1:length(maps)
    disp(['Map ' num2str(m) '...']);
    [~,V,Info,~] = BrikLoad(maps{m});
    if (nargin > 3)
        ind = find(V >= thresh(m));
        map_vals = ones(length(ind),1);
    else
        ind = find(V >= 1);
        map_vals = V(ind);
    end
    map_coords = zeros(length(ind),3);
    [map_coords(:,1), map_coords(:,2), map_coords(:,3)] = ind2sub(size(V), ind);
    origin = Info.ORIGIN;
    step_size = Info.DELTA;
    map_coords = bsxfun(@plus, map_coords*diag(step_size), origin);
    %map_coords(:,1:2) = -1*map_coords(:,1:2);
    
    [nearest,dists] = dsearchn(double(cell2mat(coords)),T,map_coords);
    nearest(dists > max_dist) = 0;
    
    overlays{m} = cell(2,1);
    for hem = 1:2
        overlays{m}{hem} = zeros(32492, 1);

        for v = 1:32492
            matching_vox = find(nearest == (v + 32492*(hem-1)));
            if (~isempty(matching_vox))
                overlays{m}{hem}(v) = mean(map_vals(matching_vox)); %mode(map_vals(matching_vox));
            end
        end
    end
end

end

