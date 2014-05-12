function M = LoadMetric(filename)

fid = fopen(filename);
xmlstr = fread(fid,'*char')';
fclose(fid);
block_starts = strfind(xmlstr, '<Data>');
block_starts = block_starts + 6;
block_ends = strfind(xmlstr, '</Data>');
block_ends = block_ends - 1;

if (length(block_starts)>1)
    disp('Multiple Data blocks');
end
M = typecast(dunzip(base64decode(xmlstr(block_starts(1):block_ends(1)))), 'single');
end