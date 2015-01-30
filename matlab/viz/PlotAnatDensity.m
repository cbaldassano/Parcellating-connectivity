function PlotAnatDensity()
load('/data/supervoxel/output/group468/totalconn_subj.mat');

figure('Color',[1 1 1],'Position',[2619         817         928         480]);

h = plot((1:20)-0.5,squeeze(mean(subj_totalconn(:,:,1:20),1))');
cmap = PTPalette(8);
cmap = cmap([4 7 5 1 6 2 3 8],:);
for i = 1:8
    set(h(i),'Color',cmap(i,:),'LineWidth',2);
end

hold on;
for i = 1:8
    patch([(1:20)-0.5 fliplr((1:20)-0.5)],...
    [squeeze(mean(subj_totalconn(:,i,1:20)))+squeeze(std(subj_totalconn(:,i,1:20)))/sqrt(10);flipud(squeeze(mean(subj_totalconn(:,i,1:20)))-squeeze(std(subj_totalconn(:,i,1:20)))/sqrt(10))],...
    cmap(i,:),'EdgeColor','none','FaceAlpha',0.4);
end

box off;
xlabel('Euclidean distance to connected voxel (cm)');
ylabel('Total structural connectivity');
legend('TOS','cIPL1','cIPL2','cIPL3','RSC','PHC1','PHC2','aPPA');

subj_id = repmat((1:10)',[1 8 20]);
cIPL_fac = repmat([0 1 1 1 0 0 0 0],[10 1 20]);
bin_fac = permute(repmat((1:20)',[1 10 8]),[2 3 1]);
stats = rm_anova2(subj_totalconn(:),subj_id(:),cIPL_fac(:),bin_fac(:),{'cIPL','dist'});
end