function consist_effect = PPASupervoxelDiff(C, coords, z, PPA_z, subregions)

consist_effect = cell(2,1);
for hem = 1:2
    if (hem == 1)
        counts = tabulate(z(PPA_z > 0 & (1:59412 <= 29696)'));
    else
        counts = tabulate(z(PPA_z > 0 & ~(1:59412 <= 29696)'));
    end
    counts = counts(:,2);
    [~,all_overlap] = sort(counts,'descend');
    overlap = all_overlap(1:subregions);

    mean_y = zeros(subregions,1);
    for i = 1:subregions
        mean_y(i) = mean(coords(z==overlap(i),2));
    end
    [~,post_ant] = sort(mean_y);
    overlap = overlap(post_ant);

    if (subregions == 3)
        consist_effect{hem} = ((C(overlap(1),:)-C(overlap(2),:)).*(C(overlap(2),:)-C(overlap(3),:)) > 0).*(C(overlap(1),:)-C(overlap(3),:));
        consist_effect{hem}(overlap(1:3)) = 0;
    elseif (subregions == 2)
        consist_effect{hem} = C(overlap(1),:)-C(overlap(2),:);
        consist_effect{hem}(overlap(1:2)) = 0;
    end
    
end
end