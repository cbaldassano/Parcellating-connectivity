function D = LoadCtyxCty(fname, FIPS)

N = size(FIPS,1);
D = zeros(N,N);

fid = fopen(fname);

l = fgetl(fid);
i = 0;
disp('Loading flow matrix...');
while ischar(l)
    if (mod(i,1000)==0)
        disp(['   ' num2str(i/275829*100) '%...']);
    end
    i = i+1;
    current = find(ismember(FIPS,l(2:6),'rows'), 1);
    last = find(ismember(FIPS,l(8:12),'rows'), 1);
    if (~isempty(current) && ~isempty(last))
        D(last,current) = str2double(l(374:380))/str2double(l(260:267));
    end
    l = fgetl(fid);
end
fclose(fid);

end