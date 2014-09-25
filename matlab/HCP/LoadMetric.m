function M = LoadMetric(filename)

fid = fopen(filename);
xmlstr = fread(fid,'*char')';
fclose(fid);
block_starts = strfind(xmlstr, '<Data>');
block_starts = block_starts + 6;
block_ends = strfind(xmlstr, '</Data>');
block_ends = block_ends - 1;

M = cell(length(block_starts),1);
for i = 1:length(M)
    M{i} = typecast(dunzip(base64decode(xmlstr(block_starts(i):block_ends(i)))), 'single');
end
end