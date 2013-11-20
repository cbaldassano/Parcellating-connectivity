function conn_diff = CalcPPAConnDiff(z, labels, coords, bold )

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

LOC_bold = mean(bold(labels==1 | labels==2,:),1);
TOS_bold = mean(bold(labels==3 | labels==4,:),1);
RSC_bold = mean(bold(labels==5 | labels==6,:),1);
IPL_bold = mean(bold(labels==7 | labels==8,:),1);
aPPA_bold = mean(bold(z==aPPA(1) | z==aPPA(2),:),1);
pPPA_bold = mean(bold(z==pPPA(1) | z==pPPA(2),:),1);

aPPA_conn = [corr2(aPPA_bold,LOC_bold) ...
             corr2(aPPA_bold,TOS_bold) ...
             corr2(aPPA_bold,RSC_bold) ...
             corr2(aPPA_bold,IPL_bold)];
pPPA_conn = [corr2(pPPA_bold,LOC_bold) ...
             corr2(pPPA_bold,TOS_bold) ...
             corr2(pPPA_bold,RSC_bold) ...
             corr2(pPPA_bold,IPL_bold)];

conn_diff = pPPA_conn - aPPA_conn;


end

