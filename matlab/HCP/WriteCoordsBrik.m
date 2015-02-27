function WriteCoordsBrik(templateFile, maxdist, vals, saveFile)

load('../../data/group468/coords');
coords = double(coords);

Optin.Frames = 1;
[err, V, Info, ErrMessage] = BrikLoad(templateFile,Optin);
V = zeros(size(V));
origin = Info.ORIGIN;
step_size = Info.DELTA;
[V_inds_x, V_inds_y, V_inds_z] = meshgrid(1:size(V,1),1:size(V,2),1:size(V,3));
V_coords = bsxfun(@plus, ([V_inds_x(:) V_inds_y(:) V_inds_z(:)]-1)*diag(step_size), origin);

T = delaunayn(coords);
[nearest, dists] = dsearchn(coords,T,V_coords);
nearest(dists>maxdist) = length(vals)+1;
vals = [vals 0];
V_z = vals(nearest);

V(sub2ind(size(V),V_inds_x(:),V_inds_y(:),V_inds_z(:))) = V_z;
Opt.Prefix = saveFile;
Info.DATASET_RANK(2) = 1;
Info.BRICK_TYPES = 0;
Info.BRICK_STATS = [0 max(vals)];
Info.BRICK_FLOAT_FACS = 0;
Info.BRICK_KEYWORDS = '';
Info.BRICK_LABS = 'Parcels';
Info.HISTORY_NOTE = '';
WriteBrik(V,Info,Opt);
end

% function WriteCoordsBrik(templateFile, coords, voxHalfSide, vals, saveFile)
% 
% [err, V, Info, ErrMessage] = BrikLoad(templateFile);
% origin = Info.ORIGIN;
% step_size = Info.DELTA;
% 
% coords = round(bsxfun(@minus, coords, origin)*diag(1./step_size));
% V = zeros(size(V));
% for i = 1:size(coords,1)
%     for x = -voxHalfSide:voxHalfSide
%         for y = -voxHalfSide:voxHalfSide
%             for z = -voxHalfSide:voxHalfSide
%                 V(coords(i,1)+x,coords(i,2)+y,coords(i,3)+z) = vals(i);
%             end
%         end
%     end
% end
% 
% Opt.Prefix = saveFile;
% WriteBrik(V,Info,Opt);
% end