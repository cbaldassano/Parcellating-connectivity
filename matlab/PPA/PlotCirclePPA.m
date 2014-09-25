function PlotCirclePPA()
load('PPA_May23_10_5_2.mat');
subj = {'100307' '103515' '103818' '111312' '117122' '118932' '119833' '120212' '125525' ...
        '128632' '130013' '137128' '138231' '142828' '143325' '144226' '149337' '156637' '159239' '161731'};
    
subjz = cell(20,1);
subjLabels = cell(20,1);
subjConn = cell(20,1);
thresh = 0.22;
for j=1%:20
    disp(subj{j});
    loaded = load([subj{j} '/PPA_scales']);
    labels = loaded.scaledLabels{2};
    coords = loaded.scaledCoords{2};
    D = atanh(corr(loaded.scaledBold{2}'));
    
    maxLLind = zeros(10,1);
    for sigInd = 1:10
        maxLL = -Inf;
        for seed = 1:5
            if (results{sigInd}{seed}(j).LL(end) > maxLL)
                maxLL = results{sigInd}{seed}(j).LL(end);
                maxLLind(sigInd) = seed;
            end
        end
        KmatMin(j,sigInd) = min([length(unique(results{sigInd}{maxLLind(sigInd)}(j).z(labels==9))) length(unique(results{sigInd}{maxLLind(sigInd)}(j).z(labels==10)))]);
    end
    chosenSig = find(KmatMin(j,:)<2,1)-1;
    
    subjz{j} = results{chosenSig}{maxLLind(chosenSig)}(j).z;
    
    %Reorder to raPPA, laPPA, rpPPA, lpPPA
    
    PPAlabels = zeros(size(subjz{j}));
    
    for side = 1:2
        hemClusts = unique(subjz{j}(labels==(side+8)));
        meanY = zeros(length(hemClusts),1);
        for cl=1:length(hemClusts)
            meanY(cl) = mean(coords(subjz{j}==hemClusts(cl),2));
        end

        [~,aPPA] = max(meanY);
        aPPA = hemClusts(aPPA);
        [~,pPPA] = min(meanY);
        pPPA = hemClusts(pPPA);
        
        PPAlabels(subjz{j}==aPPA) = side;
        PPAlabels(subjz{j}==pPPA) = side+2;
    end
    
    labels = labels+4;
    voxOrder = [];
    subjLabels{j} = [];
    for k=1:4
        subjLabels{j} = [subjLabels{j};k*ones(sum(PPAlabels==k),1)];
        voxOrder = [voxOrder;find(PPAlabels==k)'];
    end
    for k=5:12
        subjLabels{j} = [subjLabels{j};k*ones(sum(labels==k),1)];
        voxOrder = [voxOrder;find(labels==k)];
    end
    reorderD = D(voxOrder,voxOrder);
    PPAInds =1:find(subjLabels{j}<=4,1,'last');
    subjConn{j} = zeros(size(reorderD));
    for v=find(subjLabels{j}>4,1,'first'):size(reorderD,2)
        [~,highestInd] = max(reorderD(PPAInds,v));
        subjConn{j}(highestInd,v) = 1;
    end

%     otherVox = find(subjLabels{j}>4,1,'first');
%     subjConn{j} = reorderD;
%     for v=1:find(subjLabels{j}<=4,1,'last')
%         [~,highestInd] = sort(reorderD(v,otherVox:end),'descend');
%         threshed = zeros(size(reorderD,2)-otherVox+1,1);
%         threshed(highestInd(1:5)) = 1;
%         subjConn{j}(v,otherVox:end) = threshed;
%     end
%     subjConn{j} = reorderD >= thresh;
end
PlotConnCircle(subjConn, subjLabels);
end
    