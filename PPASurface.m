function PPA = PPASurface(midthickness_files, PPA_file)

[~,~,Info] = BrikLoad(PPA_file);
V_dict = cell(2,1);
origin = Info.ORIGIN;
step_size = Info.DELTA;
coords = cell(2,1);
for hem = 1:2
    fid = fopen(midthickness_files{hem});
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;
    coords{hem} = typecast(dunzip(base64decode( ...
        xmlstr(block_starts(1):block_ends(1)))), 'single');
    coords{hem} = permute(reshape(coords{hem},[3 32492]),[2 1]);
    coords{hem}(:,1:2) = -1*coords{hem}(:,1:2);
    
    V_dict{hem} = zeros(Info.DATASET_DIMENSIONS(1:3));
    for v = 1:size(coords{hem},1)
        V_ind = round((coords{hem}(v,:)-origin)./step_size) + 1;
        V_dict{hem}(V_ind(1), V_ind(2), V_ind(3)) = v;
    end
end


[~,V,~,~] = BrikLoad(PPA_file);
V(:,:,100:end) = zeros(size(V(:,:,100:end)));

PPA = cell(2,1);

ind = find(V == 4 | V == 5);
PPA{1} = zeros(32492, 1);
surf_ind = V_dict{1}(ind);
ind = ind(surf_ind > 0);
[i, j, k] = ind2sub(size(V),ind);
ind = ind(j < 125);
[i, j, k] = ind2sub(size(V),ind);
ind = ind(k <= 96 | i > 155);
surf_ind = V_dict{1}(ind);
PPA{1}(surf_ind) = 1;

ind = find(V == 4 | V == 5);
PPA{2} = zeros(32492, 1);
surf_ind = V_dict{2}(ind);
ind = ind(surf_ind > 0);
[i, j, k] = ind2sub(size(V),ind);
p1 = [107.5 130 88];
p2 = [105 136 85];
p3 = [93 122 93.5];
n = cross(p1-p2,p1-p3);
ind = ind(((n(1)*i + n(2)*j + n(3)*k) < sum(n.*p1)) & ~(i < 90 & j < 120 & k < 91));
surf_ind = V_dict{2}(ind);
PPA{2}(surf_ind) = 1;


end