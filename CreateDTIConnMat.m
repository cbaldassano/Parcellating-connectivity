function CreateDTIConnMat( subj )

path = ['/mnt/chekov_scratch/HCP/' subj '/T1w/'];
disp('Mapping voxels to surfaces...');
vox_mapping = MatchFDTCoordsToSurface([path 'gm.nii.gz'], ...
    [path 'probtrackx/1/coords_for_fdt_matrix3'],...
    {[path 'fsaverage_LR32k/' subj '.L.midthickness.32k_fs_LR.surf.gii'],...
     [path 'fsaverage_LR32k/' subj '.R.midthickness.32k_fs_LR.surf.gii']});
D = LoadTractography(vox_mapping, [path '/probtrackx/']);
save(['/data/supervoxel/data/Q3/' subj '/full.mat'], 'D', '-v7.3');
D(D > 10000) = 10000;
save(['/data/supervoxel/data/Q3/' subj '/full_thresh10000.mat'], 'D', '-v7.3');
end

