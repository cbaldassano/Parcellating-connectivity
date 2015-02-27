function WriteCifti(data_64984, fout)

template = '../../data/Conte69.MyelinAndCorrThickness.32k_fs_LR.dscalar.nii';

template = ft_read_cifti(template);
template.dscalar = data_64984;
ft_write_cifti(fout, template, 'parameter', 'dscalar');

end