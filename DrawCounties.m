function DrawCounties(borders, z)
figure; hold on;
for i = 1:length(borders) % contig'
    for j = 1:length(borders{i})
        fill(borders{i}{j}(:,1),borders{i}{j}(:,2),z(i));
    end
end
end