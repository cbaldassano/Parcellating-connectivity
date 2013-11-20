function WritePPAConnTable(z, labels, coords, bold, savename)

aPPA = [0 0]; pPPA = [0 0];
for side = 1:2
    hem_clusts = unique(z(labels==(side+8)));
    meanY = zeros(length(hem_clusts),1);
    for cl=1:length(hem_clusts)
        meanY(cl) = mean(coords(z==hem_clusts(cl),2));
    end

    [~,ant_ind] = max(meanY);
    aPPA(side) = hem_clusts(ant_ind);
    [~,pos_ind] = min(meanY);
    pPPA(side) = hem_clusts(pos_ind);
end

aPPA_bold = mean(bold(z==aPPA(1) | z==aPPA(2),:),1);
pPPA_bold = mean(bold(z==pPPA(1) | z==pPPA(2),:),1);

aPPA_numvox = sum(z==aPPA(1) | z==aPPA(2));
pPPA_numvox = sum(z==pPPA(1) | z==pPPA(2));

conn_table = zeros(4,2);
for roi = 1:4
    for v = find((labels == (roi*2 - 1)) | (labels == roi*2))'
        if (corr2(bold(v,:),pPPA_bold) > corr2(bold(v,:),aPPA_bold))
            conn_table(roi,1) = conn_table(roi,1)+1;
        else
            conn_table(roi,2) = conn_table(roi,2)+1;
        end
    end
end

roi_bold = zeros(4,size(bold,2));
for roi = 1:4
    roi_bold(roi,:) = mean(bold((labels == (roi*2 - 1)) | (labels == roi*2),:),1);
end
conn_table_inv = zeros(2,4);
for v = find(z==pPPA(1) | z==pPPA(2));
    roi_corr = zeros(4,1);
    for roi = 1:4
        roi_corr(roi) = corr2(bold(v,:), roi_bold(roi,:));
    end
    [~,max_roi] = max(roi_corr);
    conn_table_inv(1,max_roi) = conn_table_inv(1,max_roi)+1;
end
for v = find(z==aPPA(1) | z==aPPA(2));
    roi_corr = zeros(4,1);
    for roi = 1:4
        roi_corr(roi) = corr2(bold(v,:), roi_bold(roi,:));
    end
    [~,max_roi] = max(roi_corr);
    conn_table_inv(2,max_roi) = conn_table_inv(2,max_roi)+1;
end

fid = fopen(savename, 'w');
fprintf(fid,'- - 1 2 3 4 5 6\n');
fprintf(fid,'- - Posterior Anterior LOC TOS RSC cIPL\n');
fprintf(fid,'1 Posterior - - %d %d %d %d\n',conn_table_inv(1,:));
fprintf(fid,'2 Anterior - - %d %d %d %d\n',conn_table_inv(2,:));
for roi = 1:4
    fprintf(fid,'%d ', roi+2);
    switch roi
        case 1
            fprintf(fid,'LOC ');
        case 2
            fprintf(fid,'TOS ');
        case 3
            fprintf(fid,'RSC ');
        case 4
            fprintf(fid,'cIPL ');
    end
    fprintf(fid,'%d %d - - - -\n',conn_table(roi,1),conn_table(roi,2));
end
