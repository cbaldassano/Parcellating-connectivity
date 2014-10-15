function PlotAnatDensity()
load('/data/supervoxel/output/group468/degreedist.mat');

figure('Color',[1 1 1],'Position',[2619         817         928         480]);
LRdegree = (anat_density([1 2 3 4 9 11 12 13],:) + anat_density([5 6 7 8 10 14 15 16],:))/2;
h = plot((1:20)-0.5,LRdegree(:,1:20)');
box off;
xlabel('Euclidean distance to connected voxel (cm)');
ylabel('Mean structural connectivity');

cmap = PTPalette(8);
cmap = cmap([4 7 5 1 6 2 3 8],:);
for i = 1:8
    set(h(i),'Color',cmap(i,:),'LineWidth',2);
end
legend('TOS','cIPL1','cIPL2','cIPL3','RSC','PHC1','PHC2','aPPA');
end