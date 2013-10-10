function ImportGIFTI( eigenmap_files, midthickness_files, savename, savename_small)

% Inspired by Guillaume Flandin <Guillaume@artefact.tk> http://www.artefact.tk/software/matlab/gifti/
eigmaps = cell(2,1);
raw_coords = cell(2,1);
tri = cell(2,1);
orig_ind = cell(2,1);
adj_hem = cell(2,1);
for hem = 1:2
    disp(['Loading eigenmap ' num2str(hem) '...']);
    fid = fopen(eigenmap_files{hem});
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    
    disp('   Extracting data...');
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;
    
    eigmaps{hem} = zeros(32492,4500, 'single');
    for i = 1:4500
        eigmaps{hem}(:,i) = typecast(dunzip(base64decode( ...
            xmlstr(block_starts(i):block_ends(i)))), 'single');
    end
    
    disp(['Loading surface ' num2str(hem) '...']);
    fid = fopen(midthickness_files{hem});
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    
    disp('   Extracting data...');
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;
    
    raw_coords{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(1):block_ends(1)))), 'single');
    
    raw_coords{hem} = permute(reshape(raw_coords{hem},[3 32492]),[2 1]);
    
    tri{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(2):block_ends(2)))), 'int32');
    
    tri{hem} = permute(reshape(tri{hem},[3 64980]),[2 1]) + 1;
    
    % Remove invalid voxels
    orig_ind{hem} = 1:32492;
    valid_vox = any(eigmaps{hem},2);
    orig_ind{hem} = orig_ind{hem}(valid_vox);
    eigmaps{hem} = eigmaps{hem}(valid_vox,:);
    raw_coords{hem} = raw_coords{hem}(valid_vox,:);
    
    new_ind = zeros(32492, 1);
    for i = 1:length(orig_ind{hem})
        new_ind(orig_ind{hem}(i)) = i;
    end
    
    % Set up adjacency
    adj_hem{hem} = cell(size(eigmaps{hem},1),1);
    for v = 1:size(eigmaps{hem},1)
        [triangles, ~] = find(tri{hem} == orig_ind{hem}(v));
        neighbors = tri{hem}(triangles,:);
        adj_hem{hem}{v} = setdiff(new_ind(unique(neighbors(:))), [0 v]);
    end
end

disp('Creating full dataset...');
left_vox = size(eigmaps{1},1);
adj_hem{2} = cellfun(@(x) x + left_vox, adj_hem{2}, 'UniformOutput', false);

adj_list = [adj_hem{1}; adj_hem{2}];
coords = cell2mat(raw_coords);
D = atanh(corr(cell2mat(eigmaps)'));

save(savename, 'D', 'adj_list', 'coords', 'orig_ind', '-v7.3');


% Randomly remove 90% of voxels, for testing
disp('Creating small dataset...');
rng(1);
valid_vox = randperm(size(coords,1),round(size(coords,1)/10));
D = D(valid_vox,valid_vox);
coords = coords(valid_vox,:);

orig_ind{1} = orig_ind{1}(valid_vox(valid_vox <= left_vox));
orig_ind{2} = orig_ind{2}(valid_vox(valid_vox > left_vox) - left_vox);

adj_list = cell(size(coords,1),1);
kNN = 5;
for v = 1:length(adj_list)
    dists = sum((coords - repmat(coords(v,:),size(coords,1),1)).^2,2);
    [~,nearest] = sort(dists);
    adj_list{v} = nearest(2:(kNN+1));
end

save(savename_small, 'D', 'adj_list', 'coords', 'orig_ind', '-v7.3');

end
