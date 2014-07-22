function C = SupervoxelConnectivityMatrixFunc( eigmaps, z )

[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

C = zeros(length(parcels),length(parcels));
meanmaps = zeros(length(parcels),size(eigmaps,2));

%disp('Computing means...');
for i = 1:length(parcels)
    meanmaps(i,:) = mean(eigmaps(parcels{i},:),1);
end

%disp('Computing matrix...');
for i = 1:length(parcels)
    for j = (i+1):length(parcels)
        C(i,j) = CorrNoCentering(meanmaps(i,:),meanmaps(j,:));
    end
end

C = C + C';

C = atanh(C);
end

