function PlotSubjVarGeneralization()

dd216 = load('../output/group468/ddCRP2000_rng1_subj_var.mat');
dd216 = dd216.subj_var;
dd172 = load('../output/group468/ddCRP3000_rng1_subj_var.mat');
dd172 = dd172.subj_var;
dd155 = load('../output/group468/ddCRP4000_rng1_subj_var.mat');
dd155 = dd155.subj_var;
dd140 = load('../output/group468/ddCRP5000_rng1_subj_var.mat');
dd140 = dd140.subj_var;

WC216 = load('../output/group468/WC216_subj_var.mat');
WC216 = WC216.WC_subj_var;
WC172 = load('../output/group468/WC172_subj_var.mat');
WC172 = WC172.WC_subj_var;
WC155 = load('../output/group468/WC155_subj_var.mat');
WC155 = WC155.WC_subj_var;
WC140 = load('../output/group468/WC140_subj_var.mat');
WC140 = WC140.WC_subj_var;

figure('Color',[1 1 1]);
plot([140 155 172 216],[mean(dd140) mean(dd155) mean(dd172) mean(dd216)],'d',...
    [140 155 172 216],[mean(WC140) mean(WC155) mean(WC172) mean(WC216)],'d');
xlim([130 240]);
box off;

[~,p] = ttest(dd140,WC140,0.05,'right');
disp(['140 p=' num2str(p)]);
[~,p] = ttest(dd155,WC155,0.05,'right');
disp(['155 p=' num2str(p)]);
[~,p] = ttest(dd172,WC172,0.05,'right');
disp(['172 p=' num2str(p)]);
[~,p] = ttest(dd216,WC216,0.05,'right');
disp(['216 p=' num2str(p)]);
end