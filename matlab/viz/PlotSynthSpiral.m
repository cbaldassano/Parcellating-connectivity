function PlotSynthSpiral()

spiral = load('../../output/synth/spiral.mat');
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
     synth_sig,mean(spiral.DC,2));
 
xlim([0 9]);
ylim([0 1.05]);
xlabel('Noise \sigma');
ylabel('NMI with ground truth');
box off;
for i = 1:length(h)
 set(h(i),'LineWidth',2);
 set(h(i),'Color',palette(i,:));
end

subplot(1,2,2);
palette = TableauPalette(5);
figure('Color', [1 1 1], 'Position', [933         855        300         290]);
h = plot(synth_sig, mean(spiral.DC_K,2));
xlim([0 9]);
ylim([2 9.05]);
xlabel('Noise \sigma');
ylabel('Inferred number of clusters');
box off;
for i = 1:length(h)
 set(h(i),'LineWidth',2);
 set(h(i),'Color',palette(i+2,:));
end


end