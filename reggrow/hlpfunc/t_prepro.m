function [Data, Mask] = t_prepro(DATA_DS, DEMEAN, VoxVarNorm)

if ndims(DATA_DS) == 4
   [dx, dy, dz, dt] = size(DATA_DS);

    if (DEMEAN)
        %display('Calculating mean and demeaning data')
        Mean            = mean(DATA_DS,4);
        DATA_DM         = DATA_DS-repmat(Mean,[1,1,1,dt]);
    else
        DATA_DM         = DATA_DS;
    end

        %display('Calculating voxel-wise variance')
        Var             = var(DATA_DM,[],4);

        %display('Calculating mask')
        Mask            = (Var/max(Var(:)))>1e-6;

    if (VoxVarNorm) 
        %display('Normalising voxel-wise variance')
        Var(~Mask)=1;
        DATA_VN         = DATA_DM./sqrt(repmat(Var,[1 1 1 dt]));
        Var(~Mask)      = 0;                            %#ok
    else 
        DATA_VN         = DATA_DM;    
    end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Reshape data matrix
    Data                = reshape(DATA_VN,[dx*dy*dz dt]);
    Data(Mask(:)==0,:)  = [];
    Data                = Data';
elseif ndims(DATA_DS) == 3
   [dx, dy, dt] = size(DATA_DS);
    
    if (DEMEAN)
        %display('Calculating mean and demeaning data')
        Mean            = mean(DATA_DS,3);
        DATA_DM         = DATA_DS-repmat(Mean,[1,1,dt]);
    else
        DATA_DM         = DATA_DS;
    end
    
        %display('Calculating voxel-wise variance')
        Var             = var(DATA_DM,[],3);
        
        %display('Calculating mask')
        Mask            = (Var/max(Var(:)))>1e-6;
    
    if (VoxVarNorm)
        %display('Normalising voxel-wise variance')
        Var(~Mask)=1;
        DATA_VN         = DATA_DM./sqrt(repmat(Var,[1 1 dt]));
        Var(~Mask)      = 0;                            %#ok
    else 
        DATA_VN         = DATA_DM;    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Reshape data matrix
    Data                = reshape(DATA_VN,[dx*dy dt]);
    Data(Mask(:)==0,:)  = [];
    Data                = Data';
elseif ndims(DATA_DS) == 2
   [dx, dt] = size(DATA_DS);
    
    if (DEMEAN)
        %display('Calculating mean and demeaning data')
        Mean            = mean(DATA_DS,2);
        DATA_DM         = DATA_DS-repmat(Mean,[1,dt]);
    else
        DATA_DM         = DATA_DS;
    end
    
        %display('Calculating voxel-wise variance')
        Var             = var(DATA_DM,[],2);
        
        %display('Calculating mask')
        Mask            = (Var/max(Var(:)))>1e-6;
    
    if (VoxVarNorm)
        %display('Normalising voxel-wise variance')
        Var(~Mask)=1;
        DATA_VN         = DATA_DM./sqrt(repmat(Var,[1 dt]));
        Var(~Mask)      = 0;                            %#ok
    else 
        DATA_VN         = DATA_DM;    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Reshape data matrix
    Data                = reshape(DATA_VN,[dx dt]);
    Data(Mask(:)==0,:)  = [];
    Data                = Data';
    
end