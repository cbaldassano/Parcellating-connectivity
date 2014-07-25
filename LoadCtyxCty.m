function [D D_MOE Pop] = LoadCtyxCty(fname, FIPS)

N = size(FIPS,1);
D = zeros(N,N);
D_MOE = zeros(N,N);
Pop = NaN(N,1);

fid = fopen(fname,'r','n','ISO-8859-1');

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
        % http://www.census.gov/hhes/migration/files/acs/county-to-county/2007-2011/2007-2011%20Migration%20Flows%20Documentation.pdf
        if (isnan(Pop(last)))
            Pop(last) = str2double(l(260:267));
        end
        D(last,current) = str2double(l(374:380))/Pop(last);
        D_MOE(last,current) = str2double(l(382:388))/Pop(last);
    end
    l = fgetl(fid);
end
fclose(fid);

end