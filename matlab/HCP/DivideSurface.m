function DivideSurface(midthickness_files, surface_template_file, z, orig_ind, adj_list, overlay, save_prefix, overlay_suffix)
% SURF_PRE = '/data/supervoxel/data/group468/surfaces/';
%DivideSurface({[SURF_PRE 'Q1-Q6_R440.L.very_inflated.32k_fs_LR.surf.gii'],[SURF_PRE 'Q1-Q6_R440.R.very_inflated.32k_fs_LR.surf.gii']}, [SURF_PRE 'surface_template.gii'], map_z, orig_ind, adj_list, mapAtlas, '/data/supervoxel/output/group468/surf/s', 'ArcaroAtlas');

load('/data/supervoxel/data/curvature.mat');
[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));
parcel_order = cell(2,1);
orig_parcels = cell(2,1);
for i = 1:length(parcels)
    if (parcels{i}(1) <= 29696)
        orig_parcels{1} = [orig_parcels{1} {orig_ind{1}(parcels{i})}];
        parcel_order{1} = [parcel_order{1} i];
    else
        orig_parcels{2} = [orig_parcels{2} {orig_ind{2}(parcels{i} - 29696)}];
        parcel_order{2} = [parcel_order{2} i];
    end
end

% Assumes overlay is in 59412 space
rng(3);
max_overlay = max(overlay);
overlay_adj = cell(max_overlay,1);
for i = 1:max_overlay
    overlay_adj{i} = unique(overlay(setdiff(unique([adj_list{overlay == i}]),unique(find(overlay == i)))));
    overlay_adj{i} = overlay_adj{i}(overlay_adj{i}>0);
end
palette = PTPalette(12);
colors = repmat(palette,ceil(max_overlay/12),1);
colors = colors(randperm(size(colors,1),max_overlay),:);
no_conflict = false;
while ~no_conflict
    no_conflict = true;
    for i = 1:max_overlay
        if (ismember(colors(i,:),colors(overlay_adj{i},:),'rows'));
            colors(i,:) = palette(randi(12),:);
            no_conflict = false;
        end
    end
end
overlay = ConvertToOrigHem(overlay, orig_ind);

%colors = [0 0 0; 1 1 1; colors];

% max_overlay = max(overlay);
% palette = PTPalette(12);
% colors = repmat(palette,ceil(max_overlay/12),1);
% colors = colors(randperm(size(colors,1),max_overlay),:);
% colors = [0 0 0; 1 1 1; colors];
% overlay = ConvertToOrigHem(overlay, orig_ind);
% 
% max_overlay = 0.25; %0.4; 
% min_overlay = 0;% -0.08;
% colors = PTcolormap(300, [min_overlay max_overlay]);
%colors = [0 0 0; 1 1 1; PTcolormap(300, [min_overlay max_overlay])];

% cmap = PTPalette(10);
% %cmap = cmap([2 3 4 6 7 8 9 5 10 1],:);
% colors = [0 0 0; 1 1 1; cmap];

% colors = BR2Dcmap(100);
% %min1 = 12; max1 = 18;
% %min2 = 12; max2 = 18;
% min1 = 0; max1 = 700;%1000;
% min2 = 0; max2 = 700;%1000;
% overlay1 = ConvertToOrigHem(overlay(1,:), orig_ind);
% overlay2 = ConvertToOrigHem(overlay(2,:), orig_ind);

% Inspired by Guillaume Flandin <Guillaume@artefact.tk> http://www.artefact.tk/software/matlab/gifti/
fid = fopen(surface_template_file);
template_xml = fread(fid,'*char')';
fclose(fid);

fid_sizelist = fopen([save_prefix '_slist.txt'], 'w');
fid_filelist = fopen([save_prefix '_flist.txt'], 'w');
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
        [~,filename] = fileparts([save_prefix '_' num2str(hem) '_' num2str(i)]);
        fprintf(fid_filelist, '''%s'', ', filename);
        parcel_tri = tri{hem}(any(ismember(tri{hem}, orig_parcels{hem}{i}), 2), :);
        if (isempty(parcel_tri))
            continue;
        end
        
        parcel_overlay = overlay{hem}(parcel_tri);
        parcel_curv = curv{hem}(parcel_tri);
        
        %gives borders
        parcel_overlay(~ismember(parcel_tri, orig_parcels{hem}{i})) = -1;

%         parcel_overlay1 = overlay1{hem}(parcel_tri);
%         parcel_overlay2 = overlay2{hem}(parcel_tri);
%         
%         %gives borders
%         parcel_overlay1(~ismember(parcel_tri, orig_parcels{hem}{i})) = -1;
        
        fid = fopen([save_prefix '_' num2str(hem) '_' num2str(i) '_' overlay_suffix], 'w');
        %for f = 1:size(parcel_overlay1,1)
        for f = 1:size(parcel_overlay,1)
            for v = 1:3
                
                val = parcel_overlay(f,v);
                if (val == -1)
                    fprintf(fid, '0 0 0\n');
                elseif (val==0)
                     if (parcel_curv(f,v)>=0)
                         fprintf(fid, '0.698 0.698 0.698\n');
                     else
                         fprintf(fid, '0.4 0.4 0.4\n');
                     end
                else
                    fprintf(fid, '%f %f %f\n', colors(val,1), colors(val,2), colors(val,3));
                end
                
                
%                 if (~ismember(parcel_tri(f,v),orig_parcels{hem}{i}))
%                    fprintf(fid, '0.000000 0.000000 0.000000\n');
%                 elseif (parcel_overlay(f,v)==-1)
%                    fprintf(fid, '0.500000 0.500000 0.500000\n');
%                 else
%                     ind = max(min(round((parcel_overlay(f,v)-min_overlay)/(max_overlay-min_overlay)*size(colors,1)) + 1, size(colors,1)),1);
%                     fprintf(fid, '%f %f %f\n', colors(ind,1), colors(ind,2), colors(ind,3));
%                 end

%                  val = parcel_overlay(f,v) + 2;
%                  fprintf(fid, '%f %f %f\n', colors(val,1), colors(val,2), colors(val,3));
                
%                if (parcel_overlay(f,v)==-2)
%                    fprintf(fid, '0 0 0\n');
%                elseif (parcel_overlay(f,v)==-1)
%                    fprintf(fid, '0.5 0.5 0.5\n');
%                else
%                    ind = max(min(round((parcel_overlay(f,v)-min_overlay)/(max_overlay-min_overlay)*300) + 1, 300),1);
%                    fprintf(fid, '%f %f %f\n', colors(ind,1), colors(ind,2), colors(ind,3));
%                end

%                if (parcel_overlay1(f,v)==-1)
%                    fprintf(fid, '0 0 0\n');
%                else
%                    ind1 = max(min(round((parcel_overlay1(f,v)-min1)/(max1-min1)*100) + 1, 100),1);
%                    ind2 = max(min(round((parcel_overlay2(f,v)-min2)/(max2-min2)*100) + 1, 100),1);
%                    fprintf(fid, '%f %f %f\n', colors(ind1,ind2,1), colors(ind1,ind2,2), colors(ind1,ind2,3));
%                end

%                if (parcel_overlay(f,v) == 1)
%                    fprintf(fid, '0.000000 0.000000 0.000000\n');
%                else
%                    fprintf(fid, '1.000000 1.000000 1.000000\n');
%               end
               
%                if (~all(parcel_overlay(f,:) == parcel_overlay(f,v)))
%                    fprintf(fid, '0.000000 0.000000 0.000000\n');
%                else
%                    fprintf(fid, '1.000000 1.000000 1.000000\n');
%                end
            end
        end
        fclose(fid);
        fprintf(fid_sizelist, '%d, ', 9*size(parcel_overlay,1));
        %fprintf(fid_sizelist, '%d, ', 9*size(parcel_overlay1,1));
        
        continue;
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

fclose(fid_sizelist);
fclose(fid_filelist);
end

