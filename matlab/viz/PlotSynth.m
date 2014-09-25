function PlotSynth()

square = load('../output/synth/square.mat');
face = load('../output/synth/face.mat');
stripes = load('../output/synth/stripes.mat');
synth_sig = linspace(0,9,10);

palette = PTPalette(12);
figure('Color', [1 1 1], 'Position', [933         855        1060         290]);
subplot(1,3,1);
patch([synth_sig fliplr(synth_sig)],...
 [mean(square.RC,2)+std(square.RC,[],2);flipud(mean(square.RC,2)-std(square.RC,[],2))],...
 [0.7 0.7 0.7],'EdgeColor','none');
hold on;
h = plot(synth_sig,mean(square.LS,2),...
     synth_sig,mean(square.NC,2),...
     synth_sig,mean(square.RG,2),...
     synth_sig,mean(square.WC,2),...
     synth_sig,mean(square.DC,2));
 
xlim([0 9]);
ylim([0 1.05]);
xlabel('Noise \sigma');
ylabel('NMI with ground truth');
box off;
for i = 1:length(h)
 set(h(i),'LineWidth',2);
 set(h(i),'Color',palette(i,:));
end


subplot(1,3,2);
patch([synth_sig fliplr(synth_sig)],...
[mean(face.RC,2)+std(face.RC,[],2);flipud(mean(face.RC,2)-std(face.RC,[],2))],...
[0.7 0.7 0.7],'EdgeColor','none');
hold on;
h = plot(synth_sig,mean(face.LS,2),...
 synth_sig,mean(face.NC,2),...
 synth_sig,mean(face.RG,2),...
 synth_sig,mean(face.WC,2),...
 synth_sig,mean(face.DC,2));

xlim([0 9]);
ylim([0 1.05]);
xlabel('Noise \sigma');
ylabel('NMI with ground truth');
box off;
for i = 1:length(h)
 set(h(i),'LineWidth',2);
 set(h(i),'Color',palette(i,:));
end

subplot(1,3,3);
patch([synth_sig fliplr(synth_sig)],...
[mean(stripes.RC,2)+std(stripes.RC,[],2);flipud(mean(stripes.RC,2)-std(stripes.RC,[],2))],...
[0.7 0.7 0.7],'EdgeColor','none');
hold on;
h = plot(synth_sig,mean(stripes.LS,2),...
 synth_sig,mean(stripes.NC,2),...
 synth_sig,mean(stripes.RG,2),...
 synth_sig,mean(stripes.WC,2),...
 synth_sig,mean(stripes.DC,2));

xlim([0 9]);
ylim([0 1.05]);
xlabel('Noise \sigma');
ylabel('NMI with ground truth');
box off;
for i = 1:length(h)
 set(h(i),'LineWidth',2);
 set(h(i),'Color',palette(i,:));
end


palette = TableauPalette(5);
figure('Color', [1 1 1], 'Position', [933         855        300         290]);
h = plot(synth_sig, mean(square.DC_K,2),...
    synth_sig, mean(face.DC_K,2),...
    synth_sig, mean(stripes.DC_K,2));
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