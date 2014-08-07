function PSCO = LoadPVS(fname, Pop_orig, FIPS)
raw = importdata(fname);
raw.textdata = raw.textdata(2:end,:);

PSCO = zeros(length(FIPS),4);
for i = 1:length(FIPS)
    if (mod(i,100)==1)
        disp(num2str(i));
    end
    totalPop = 0;
    for j = 1:length(FIPS{i})
        ind = find(cellfun(@(x) strcmp(x,FIPS{i}{j}), raw.textdata(:,6)),1);
        if (isempty(ind))
            disp(['No match for ' FIPS{i}{j} ', at i=' num2str(i) ', j=' num2str(j)]);
        else
            PSCO(i,:) = PSCO(i,:) + Pop_orig{i}{j}*raw.data(ind,7:10);
            totalPop = totalPop + Pop_orig{i}{j};
        end
    end
    PSCO(i,:) = PSCO(i,:)/totalPop;
end
end