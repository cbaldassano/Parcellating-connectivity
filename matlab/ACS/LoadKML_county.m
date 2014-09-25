function [FIPS names borders adj_list] = LoadKML_county(fname)

disp('Extracting data...');
fid = fopen(fname);
xmlstr = fread(fid,'*char')';
fclose(fid);

statefp = strfind(xmlstr, ['<th>STATEFP</th>' char(13) char(10) '<td>']) + ...
    length(['<th>STATEFP</th>' char(13) char(10) '<td>']);
countyfp = strfind(xmlstr, ['<th>COUNTYFP</th>' char(13) char(10) '<td>']) + ...
    length(['<th>COUNTYFP</th>' char(13) char(10) '<td>']);

N = length(statefp);

FIPS = char(zeros(N,5));
for i = 1:N
    FIPS(i,:) = [xmlstr(statefp(i):(statefp(i)+1)) xmlstr(countyfp(i):(countyfp(i)+2))];
end

names = regexp(xmlstr, ['<name>[^<]*</name>' char(13) char(10)],'match');
names = names(3:end); %Remove header info
for i = 1:N
    names{i} = names{i}(7:(end-9));
end

coord_ind = strfind(xmlstr, '<coordinates>') + length('<coordinates>');
coord_ind_end = strfind(xmlstr, '</coordinates>')-1;
coordinates = cell(length(coord_ind),1);
for i = 1:length(coord_ind)
    collected = textscan(xmlstr(coord_ind(i):coord_ind_end(i)), '%f,%f,0','CollectOutput',1);
    coordinates{i} = collected{1};
end

borders = cell(N,1);
entry_lim = [statefp Inf];
for i = 1:N
    borders{i} = coordinates((coord_ind > entry_lim(i)) & (coord_ind < entry_lim(i+1)));
end

disp('Computing adjacency...');
borders_collapsed = cellfun(@(x) cell2mat(x), borders, 'UniformOutput', false);
adj_list = cell(N,1);
parfor i = 1:N
    adj_list{i} = setdiff(find(cellfun(@(x) any(ismember(x, borders_collapsed{i}, 'rows')), borders_collapsed)),i);
end
end