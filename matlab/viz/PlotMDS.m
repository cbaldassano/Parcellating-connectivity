function PlotMDS()

subj = load('/data/supervoxel/output/group468/ddCRP3000_dist_mds_subjprocrustes');
group = load('/data/supervoxel/output/group468/ddCRP3000_dist_mds');

paths = [...
    55    54    60    61; ...
    141   140   145   146; ...
    8    17    16    0; ...
    93   101   100 0];
RSC = [13 98];

figure('Color',[1 1 1],'Position',[676         340        1393/2        1101/2]);
scatter3(group.Y(:,1),group.Y(:,2),group.Y(:,3),20,group.Ycol,'filled'); hold on;
PlotPaths(paths, RSC, group.Y,group.Ycol,3);
scatter3(group.Y(paths(paths>0),1),group.Y(paths(paths>0),2),group.Y(paths(paths>0),3),50,group.Ycol(paths(paths>0),:),'filled');
view(-20,31);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ZTickLabel',[]);
rgb2cm;

figure('Color',[1 1 1],'Position',[676         340+1101/2        1393/2        1101/2]);
for s = 1:20
    PlotPaths(paths, RSC, subj.Ysubj{s},group.Ycol,2);
    view(-20,31);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ZTickLabel',[]);
end
rgb2cm;

figure('Color',[1 1 1],'Position',[676+1393/2         340        1393/2        1101/2]);
dorsalConnLeft = zeros(4,20);
dorsalConnRight = zeros(4,20);
for s = 1:20
    dorsalConnLeft(:,s) = subj.Csubj{s}(RSC(1),paths(1,:));
    dorsalConnRight(:,s) = subj.Csubj{s}(RSC(2),paths(2,:));
end

disp('Dorsal: ');
disp('   Left: ');
for i = 1:3
    t = tstat(dorsalConnLeft(i,:)-dorsalConnLeft(i+1,:));
    [~,p] = ttest(dorsalConnLeft(i,:),dorsalConnLeft(i+1,:));
    disp(['      ' num2str(i) ' vs ' num2str(i+1) ': ' num2str(p) ', ' num2str(t)]);
end
disp('   Right: ');
for i = 1:3
    t = tstat(dorsalConnRight(i,:)-dorsalConnRight(i+1,:));
    [~,p] = ttest(dorsalConnRight(i,:),dorsalConnRight(i+1,:));
    disp(['      ' num2str(i) ' vs ' num2str(i+1) ': ' num2str(p) ', ' num2str(t)]);
end
linecol = PTPalette(2);
%plot(1:4,dorsalConnLeft,'Color',linecol(1,:),'LineWidth',1); hold on;
%plot((1:4)+0.3,dorsalConnRight,'Color',linecol(1,:),'LineWidth',1);
plot(1:4,mean(dorsalConnLeft,2),'--o','Color',linecol(1,:),'LineWidth',3); hold on;
plot((1:4)+0.3,mean(dorsalConnRight,2),'--o','Color',linecol(1,:),'LineWidth',3);
plot(1:4,group.C(RSC(1),paths(1,:)),'-o','Color',linecol(1,:),'LineWidth',3); hold on;
plot((1:4)+0.3,group.C(RSC(2),paths(2,:)),'-o','Color',linecol(1,:),'LineWidth',3);
ylabel('Connectivity with RSC parcel');
ylim([0.2 0.75]);
set(gca,'XTick',[]);
box off;

figure('Color',[1 1 1],'Position',[676+1393/2         340+1101/2        1393/2        1101/2]);
ventralConnLeft = zeros(3,20);
ventralConnRight = zeros(3,20);
for s = 1:20
    ventralConnLeft(:,s) = subj.Csubj{s}(paths(1,4),paths(3,1:3));
    ventralConnRight(:,s) = subj.Csubj{s}(paths(2,4),paths(4,1:3));
end
disp('Ventral: ');
disp('   Left: ');
for i = 1:2
     t = tstat(ventralConnLeft(i,:)-ventralConnLeft(i+1,:));
    [~,p] = ttest(ventralConnLeft(i,:),ventralConnLeft(i+1,:));
    disp(['      ' num2str(i) ' vs ' num2str(i+1) ': ' num2str(p) ', ' num2str(t)]);
end
disp('   Right: ');
for i = 1:2
     t = tstat(ventralConnRight(i,:)-ventralConnRight(i+1,:));
    [~,p] = ttest(ventralConnRight(i,:),ventralConnRight(i+1,:));
    disp(['      ' num2str(i) ' vs ' num2str(i+1) ': ' num2str(p) ', ' num2str(t)]);
end
linecol = PTPalette(2);
%plot(1:3,ventralConnLeft,'Color',linecol(2,:),'LineWidth',1); hold on;
%plot((1:3)+0.3,ventralConnRight,'Color',linecol(2,:),'LineWidth',1);
plot(1:3,mean(ventralConnLeft,2),'--o','Color',linecol(2,:),'LineWidth',3); hold on;
plot((1:3)+0.3,mean(ventralConnRight,2),'--o','Color',linecol(2,:),'LineWidth',3);
plot(1:3,group.C(paths(1,4),paths(3,1:3)),'-o','Color',linecol(2,:),'LineWidth',3); hold on;
plot((1:3)+0.3,group.C(paths(2,4),paths(4,1:3)),'-o','Color',linecol(2,:),'LineWidth',3);
ylabel('Connectivity with cIPL3 parcel');
ylim([0.2 0.75]);
set(gca,'XTick',[]);
box off;
end

function PlotPaths(paths, RSC, Y, Ycol, width)

scatter3(Y(RSC,1),Y(RSC,2),Y(RSC,3),50,Ycol(RSC,:),'filled');
hold on;

linecol = PTPalette(2);
for i = 1:3
    line([Y(paths(1,i),1) Y(paths(1,i+1),1)],...
         [Y(paths(1,i),2) Y(paths(1,i+1),2)],...
         [Y(paths(1,i),3) Y(paths(1,i+1),3)],...
         'Color',linecol(1,:),'LineWidth',width);
end
for i = 1:3
    line([Y(paths(2,i),1) Y(paths(2,i+1),1)],...
         [Y(paths(2,i),2) Y(paths(2,i+1),2)],...
         [Y(paths(2,i),3) Y(paths(2,i+1),3)],...
         'Color',linecol(1,:),'LineWidth',width);
end
for i = 1:2
    line([Y(paths(3,i),1) Y(paths(3,i+1),1)],...
         [Y(paths(3,i),2) Y(paths(3,i+1),2)],...
         [Y(paths(3,i),3) Y(paths(3,i+1),3)],...
         'Color',linecol(2,:),'LineWidth',width);
end
for i = 1:2
    line([Y(paths(4,i),1) Y(paths(4,i+1),1)],...
         [Y(paths(4,i),2) Y(paths(4,i+1),2)],...
         [Y(paths(4,i),3) Y(paths(4,i+1),3)],...
         'Color',linecol(2,:),'LineWidth',width);
end


end

function t = tstat(x)
t = mean(x)/(std(x)/sqrt(length(x)));
end