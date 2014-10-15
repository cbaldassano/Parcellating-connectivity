function [vertices, sources, meta_coords] = MetaToVertex(filename, coords)

coords(:,1:2) = -1*coords(:,1:2);

T = delaunayn(double(coords));

fid = fopen(filename, 'r');
data = textscan(fid, '%s\t%d,%d,%d\t%s\n');
fclose(fid);
sources = data{1};

meta_coords = double([data{2} data{3} data{4}]);
tal_coords = cellfun(@(x) x == 'T', data{5});
meta_coords(tal_coords,:) = tal2mni(meta_coords(tal_coords,:));
[vertices,dists] = dsearchn(double(coords),T,meta_coords);

end