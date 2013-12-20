function ClusteringMatrix(D, parcels)

v_ordered = [];
for i = 1:length(parcels);
    v_ordered = [v_ordered parcels{i}(randperm(length(parcels{i})))];
end
v_unordered = v_ordered(randperm(length(v_ordered)));

figure('Position',[100 1000 560 560], 'Color', [1 1 1]);
imagesc(D(v_unordered, v_unordered));
caxis([0 1]);
axis square;
set(gca,'XTick',[]);
set(gca,'YTick',[]);

figure('Position',[700 1000 560 560], 'Color', [1 1 1]);
imagesc(D(v_ordered, v_ordered));
caxis([0 1]);
axis square;
set(gca,'XTick',[]);
set(gca,'YTick',[]);

end