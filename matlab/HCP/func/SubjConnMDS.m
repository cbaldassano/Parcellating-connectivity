function [Ysubj Csubj] = SubjConnMDS( z, groupY)
subj = {'120515','113922','111413','109325','109123','108828','108525','108323','108121','107422','107321','106521','105014','104820','103111','102311','102008','101309','101107','101006'};

[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

Ysubj = cell(length(subj),1);
Csubj = cell(length(subj),1);
for s = 1:length(subj)
    disp(['Subj ' subj{s} '...']);
    load(['/data/supervoxel/data/S500/' subj{s} '/bold']);
    Csubj{s} = zeros(length(parcels),length(parcels));
    meanmaps = zeros(length(parcels),size(bold,2));
    
    for i = 1:length(parcels)
        meanmaps(i,:) = mean(bold(parcels{i},:),1);
    end

    for i = 1:length(parcels)
        for j = (i+1):length(parcels)
            Csubj{s}(i,j) = corr2(meanmaps(i,:),meanmaps(j,:));
        end
    end

    Csubj{s} = Csubj{s} + Csubj{s}';
    Csubj{s} = atanh(Csubj{s});
    
    dist = max(Csubj{s}(:)) - Csubj{s};
    dist(1:(size(dist,1)+1):end) = 0;
    dist = squareform(dist);
    Ysubj{s} = cmdscale(dist);
    [~,Ysubj{s}] = procrustes(groupY(:,1:3), Ysubj{s}(:,1:3));
end
    
end

