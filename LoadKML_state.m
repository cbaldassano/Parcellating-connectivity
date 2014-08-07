function [FIPS borders] = LoadKML_state(fname, z_state)

disp('Extracting data...');
fid = fopen(fname);
xmlstr = fread(fid,'*char')';
fclose(fid);

statefp = strfind(xmlstr, ['<th>STATEFP</th>' char(13) char(10) '<td>']) + ...
    length(['<th>STATEFP</th>' char(13) char(10) '<td>']);

N = length(statefp);

FIPS = char(zeros(N,2));
for i = 1:N
    FIPS(i,:) = xmlstr(statefp(i):(statefp(i)+1));
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

FIPS_num = zeros(N,1);
for i = 1:N
    FIPS_num(i) = str2double(FIPS(i,:));
end
FIPS = FIPS(ismember(FIPS_num,unique(z_state)),:);
borders = borders(ismember(FIPS_num,unique(z_state)));
end