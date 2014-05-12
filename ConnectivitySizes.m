function [C_func, C_anat_mean, C_anat_alpha1, C_anat_max] = ConnectivitySizes(eigmaps, D_anat, Z, sizes)

C_func = cell(length(sizes),1);
C_anat_mean = cell(length(sizes),1);
C_anat_alpha1 = cell(length(sizes),1);
C_anat_max = cell(length(sizes),1);
for i = 1:length(sizes)
    disp(num2str(i));
    z = cluster(Z, 'maxclust', sizes(i))';
    C_func{i} = SupervoxelConnectivityMatrixFunc(eigmaps, z);
    C_anat_mean{i} = SupervoxelConnectivityMatrixAnat(D_anat, z);
    C_anat_alpha1{i} = SupervoxelConnectivityMatrixAnat(D_anat, z, 1);
    C_anat_max{i} = SupervoxelConnectivityMatrixAnat(D_anat, z, Inf);
end
end