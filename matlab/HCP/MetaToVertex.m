function meta = MetaToVertex(filename, coords)

fid = fopen(filename, 'r');
while(~feof(fid))
    category = fgetl(fid);
    category(ismember(category,' -')) = [];
    data = textscan(fid, '%s\t%d,%d,%d\t%s\n');
    num_entries = length(data{1})-1;
    meta.(category).source = data{1}(1:num_entries);
    meta.(category).metacoords = double([data{2}(1:num_entries) data{3}(1:num_entries) data{4}(1:num_entries)]);
    meta.(category).coordtype = data{5}(1:num_entries);
    fgetl(fid);
end
fclose(fid);

coords(:,1:2) = -1*coords(:,1:2);
coords = double(coords);
T = delaunayn(coords);

fields = fieldnames(meta);
for i = 1:numel(fields)
    tal_coords = cellfun(@(x) x == 'T', meta.(fields{i}).coordtype);
    meta.(fields{i}).metacoords(tal_coords,:) = tal2mni(meta.(fields{i}).metacoords(tal_coords,:));
    [meta.(fields{i}).vertices,dists] = dsearchn(coords,T,meta.(fields{i}).metacoords);
end




end