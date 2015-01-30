function PlotSynthSpiral()

spiral = load('../../output/synth/spiral.mat');
DC_originit = load('../../output/synth/spiral100.mat');
synth_sig = linspace(0,9,10);

palette = PTPalette(12);
figure('Color', [1 1 1], 'Position', [933         855        1060         290]);
subplot(1,2,1);
patch([synth_sig fliplr(synth_sig)],...
 [mean(spiral.RC,2)+std(spiral.RC,[],2);flipud(mean(spiral.RC,2)-std(spiral.RC,[],2))],...
 [0.7 0.7 0.7],'EdgeColor','none');
hold on;
h = plot(synth_sig,mean(spiral.LS,2),...
     synth_sig,mean(spiral.NC,2),...
     synth_sig,mean(spiral.RG,2),...
     synth_sig,mean(spiral.WC,2),...
     synth_sig,mean(DC_originit.DC,2));
 
xlim([0 9]);
ylim([0 1.05]);
xlabel('Noise \sigma');
ylabel('NMI with ground truth');
box off;
for i = 1:length(h)
 set(h(i),'LineWidth',2);
 set(h(i),'Color',palette(i,:));
end


DC_randinit = load('../../output/synth/spiral_randinit1000');
DC_ward10 = load('../../output/synth/spiral_ward10_100');
DC_ward2 = load('../../output/synth/spiral_ward2_100');

subplot(1,2,2);
palette = TableauPalette(5);
hold on;
 
patch([synth_sig fliplr(synth_sig)],...
 [mean(DC_originit.DC,2)+std(DC_originit.DC,[],2);flipud(mean(DC_originit.DC,2)-std(DC_originit.DC,[],2))],...
 [0 0 0],'EdgeColor','none','FaceAlpha',0.4);

patch([synth_sig fliplr(synth_sig)],...
 [mean(DC_randinit.DC,2)+std(DC_randinit.DC,[],2);flipud(mean(DC_randinit.DC,2)-std(DC_randinit.DC,[],2))],...
 palette(2,:),'EdgeColor','none','FaceAlpha',0.4);

patch([synth_sig fliplr(synth_sig)],...
 [mean(DC_ward10.DC,2)+std(DC_ward10.DC,[],2);flipud(mean(DC_ward10.DC,2)-std(DC_ward10.DC,[],2))],...
 palette(3,:),'EdgeColor','none','FaceAlpha',0.4);

patch([synth_sig fliplr(synth_sig)],...
 [mean(DC_ward2.DC,2)+std(DC_ward2.DC,[],2);flipud(mean(DC_ward2.DC,2)-std(DC_ward2.DC,[],2))],...
 palette(5,:),'EdgeColor','none','FaceAlpha',0.4);

h = plot(synth_sig,mean(DC_originit.DC,2),'--',...
     synth_sig,mean(DC_randinit.DC,2),...
     synth_sig,mean(DC_ward10.DC,2),...
     synth_sig,mean(DC_ward2.DC,2));
set(h(1),'LineWidth',2);
set(h(1),'Color',[0 0 0]);
set(h(2),'LineWidth',2);
set(h(2),'Color',palette(2,:));
set(h(3),'LineWidth',2);
set(h(3),'Color',palette(3,:));
set(h(4),'LineWidth',2);
set(h(4),'Color',palette(5,:));
 
xlim([0 9]);
ylim([0 1.05]);
xlabel('Noise \sigma');
ylabel('NMI with ground truth');
box off;


end