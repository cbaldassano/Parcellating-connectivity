function DrawCounties(borders, contig, z)
figure; hold on;
for i = contig'
    for j = 1:length(borders{i})
        fill(borders{i}{j}(:,1),borders{i}{j}(:,2),z(i));
    end
end
end