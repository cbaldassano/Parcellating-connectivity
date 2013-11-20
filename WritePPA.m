function WritePPA( subject, pos, z )

    [err,Info] = BrikInfo(['/connectome/chrisb33/' subject '/rest1_LR_bp+tlrc.BRIK']);
    Info.DATASET_RANK(2) = 1;
    Info.TAXIS_NUMS = [];
    Info.TAXIS_FLOATS = [];
    Info.TAXIS_OFFSETS = [];
    Info.BRICK_TYPES = 3;
    Info.BRICK_FLOAT_FACS = 0;
    Info.BRICK_LABS = 'clusterID';
    Info.TYPESTRING = '3DIM_HEAD_FUNC';
    Info.SCENE_DATA = [2 11 1 -999 -999 -999 -999 -999];

    Xmax = Info.DATASET_DIMENSIONS(1);
    Ymax = Info.DATASET_DIMENSIONS(2);
    Zmax = Info.DATASET_DIMENSIONS(3);
    M = zeros(Xmax,Ymax,Zmax);
    
    for i = 1:size(pos,1)
        M(pos(i,1)+1,pos(i,2)+1,pos(i,3)+1)= z(i);
    end

    Info.BRICK_STATS = [min(z) max(z)];
    
    Opt.Prefix = 'PPAclusts';
    delete([Opt.Prefix '*']);
    [err, ErrMessage, Info] = WriteBrik(M,Info,Opt);


end

