function [w_col overlay] = SubjSvmWeights(models, parcel_inds)

subj_w = -1*cell2mat(cellfun(@(x) x.SVs'*x.sv_coef, models, 'UniformOutput', false)');
mean_w = mean(subj_w,2);
[~,p] = ttest(subj_w',0,0.05,'right');
sig = p<0.05;

minw = 0;%-0.4;
maxw = 0.4;

numc = 100;
cmap = PTcolormap(numc,[minw maxw]);
col_ind =  round((mean_w-minw)/(maxw-minw)*numc);
col_ind(col_ind<1) = 1;
col_ind(col_ind>100) = 100;
w_col = cmap(col_ind,:);
w_col(~sig,:) = repmat([0.5 0.5 0.5],sum(~sig),1);

overlay = mean_w;
overlay(~sig) = -1;

for i = parcel_inds
    [~,p] = ttest(subj_w(i,:),0,0.05,'right');
    disp(['Parcel ' num2str(i) ': t=' num2str(tstat(subj_w(i,:))) ', p=' num2str(p)]);
end
end

function t = tstat(x)
t = mean(x)/(std(x)/sqrt(length(x)));
end