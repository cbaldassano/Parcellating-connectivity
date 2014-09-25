figure('Color', [1 1 1]);

subplot(2,3,1);
smoothhist2D([sum_func sum_logD]./59411,10,[500 500],0);
title(['r = ' num2str(corr2(sum_func, sum_logD))]);
xlabel('Mean functional conn');
ylabel('Mean anatomical conn');
set(gca,'YDir','Normal');
cRange = get(gca,'clim');
colormap(PTcolormap(cRange(2),cRange));
%xlim([0 0.25]);
ylim([0 4]);

subplot(2,3,2);
smoothhist2D([sum_func./59411 curv],10,[500 500],0);
title(['r = ' num2str(corr2(sum_func, curv))]);
xlabel('Mean functional conn');
ylabel('Curvature');
set(gca,'YDir','Normal');
cRange = get(gca,'clim');
colormap(PTcolormap(cRange(2),cRange));
%xlim([0 0.25]);
ylim([-0.3 0.3]);

subplot(2,3,3);
smoothhist2D([sum_logD./59411 curv],10,[500 500],0);
title(['r = ' num2str(corr2(sum_logD, curv))]);
xlabel('Mean anatomical conn');
ylabel('Curvature');
set(gca,'YDir','Normal');
cRange = get(gca,'clim');
colormap(PTcolormap(cRange(2),cRange));
xlim([0 4]);
ylim([-0.3 0.3]);

subplot(2,3,4);
smoothhist2D([corr_map(~isnan(corr_map)),curv(~isnan(corr_map))],10,[500 500],0);
title(['r = ' num2str(corr2(corr_map(~isnan(corr_map)),curv(~isnan(corr_map))))]);
xlabel('Functional/Anatomical Correlation');
ylabel('Curvature');
set(gca,'YDir','Normal');
cRange = get(gca,'clim');
colormap(PTcolormap(cRange(2),cRange));
xlim([0 0.5]);
ylim([-0.3 0.3]);

subplot(2,3,5);
smoothhist2D([corr_map(~isnan(corr_map)),sum_func(~isnan(corr_map))./59411],10,[500 500],0);
title(['r = ' num2str(corr2(corr_map(~isnan(corr_map)),sum_func(~isnan(corr_map))))]);
xlabel('Functional/Anatomical Correlation');
ylabel('Mean functional conn');
set(gca,'YDir','Normal');
cRange = get(gca,'clim');
colormap(PTcolormap(cRange(2),cRange));
xlim([0 0.5]);
%ylim([0 0.25]);

subplot(2,3,6);
smoothhist2D([corr_map(~isnan(corr_map)),sum_logD(~isnan(corr_map))./59411],10,[500 500],0);
title(['r = ' num2str(corr2(corr_map(~isnan(corr_map)),sum_logD(~isnan(corr_map))))]);
xlabel('Functional/Anatomical Correlation');
ylabel('Mean anatomical conn');
set(gca,'YDir','Normal');
cRange = get(gca,'clim');
colormap(PTcolormap(cRange(2),cRange));
xlim([0 0.5]);
ylim([0 4]);