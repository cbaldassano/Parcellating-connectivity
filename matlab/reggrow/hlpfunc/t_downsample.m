function[DATA_DS] = t_downsample(DATA_o, DOWNSAMPLE)

if ndims(DATA_o) == 4
    if (DOWNSAMPLE~=0)
        ds = DOWNSAMPLE;
        display(['Downsampling data by averaging over patches of ' num2str(ds)  ' x ' num2str(ds) ' x ' num2str(ds)  ' voxels...'])
        [dx, dy, dz, dt] = size(DATA_o);
        DATA  = zeros(ds*ceil(dx/ds),ds*ceil(dy/ds),ds*ceil(dz/ds),dt);
        DATA(1:dx,1:dy,1:dz,:) = DATA_o;
        dx=ceil(dx/ds);
        dy=ceil(dy/ds);
        dz=ceil(dz/ds);
        DATA_DS = zeros(dx,dy,dz,dt);
        for ix = 1 : dx
            % figure(1);subplot(1,2,1);imagesc(squeeze(DATA(round(ds*ix-ds/2),:,:,50)))
            for iy = 1 : dy
                for iz = 1 : dz
                    M                   = DATA(1+(ix-1)*ds:ix*ds, 1+(iy-1)*ds:iy*ds, 1+(iz-1)*ds:iz*ds, :);
                    rM                  = reshape(M,[ds*ds*ds dt]);
                    DATA_DS(ix,iy,iz,:) = reshape(mean(rM,1),[1 1 1 dt]);
                end
            end
            % figure(1);subplot(1,2,2);imagesc(squeeze(DATA_DS(ix,:,:,50))),title([num2str(ix) 'x' num2str(iy) 'x' num2str(iz)])
            % pause(0.1)
        end
        clear DATA;
        clear ds;
        display('done')
    else
        DATA_DS = DATA_o;
        %[dx, dy, dz, dt] = size(DATA_o);
    end
elseif ndims(DATA_o) == 3
    if (DOWNSAMPLE~=0)
        ds = DOWNSAMPLE;
        display(['Downsampling data by averaging over patches of ' num2str(ds)  ' x ' num2str(ds)  ' pixels...'])
        [dx, dy, dt] = size(DATA_o);
        DATA  = zeros(ds*ceil(dx/ds),ds*ceil(dy/ds),dt);
        DATA(1:dx,1:dy,:) = DATA_o;
        dx=ceil(dx/ds);
        dy=ceil(dy/ds);
        DATA_DS = zeros(dx,dy,dt);
        for ix = 1 : dx
            % figure(1);subplot(1,2,1);imagesc(squeeze(DATA(round(ds*ix-ds/2),:,:,50)))
            for iy = 1 : dy
                M                   = DATA(1+(ix-1)*ds:ix*ds, 1+(iy-1)*ds:iy*ds, :);
                rM                  = reshape(M,[ds*ds dt]);
                DATA_DS(ix,iy,iz,:) = reshape(mean(rM,1),[1 1 dt]);
            end
            % figure(1);subplot(1,2,2);imagesc(squeeze(DATA_DS(ix,:,:,50))),title([num2str(ix) 'x' num2str(iy) 'x' num2str(iz)])
            % pause(0.1)
        end
        clear DATA;
        clear ds;
        display('done')
    else
        DATA_DS = DATA_o;
    end
elseif ndims(DATA_o)==2
    if (DOWNSAMPLE~=0)
        ds = DOWNSAMPLE;
        display(['Downsampling data by averaging over patches of ' num2str(ds) 'elements...'])
        [dx, dt] = size(DATA_o);
        DATA = zeros(ds*ceil(dx/ds),dt);
        DATA(1:dx,:) = DATA_o;
        dx=ceil(dx/ds);
        DATA_DS = zeros(dx,dt);
        for ix = 1 : dx
%                     M             = DATA(1+(ix-1)*ds:ix*ds, 1+(iy-1)*ds:iy*ds, 1+(iz-1)*ds:iz*ds, :);
%                     rM            = reshape(M,[ds*ds*ds dt]);
%                     DATA_DS(ix,:) = reshape(mean(rM,1),[1 dt]);
            DATA_DS(ix,:) = mean(DATA(1+(ix-1)*ds:ix*ds, :),1);
        end
        clear DATA;
        clear ds;
        display('done')
        
    else
        DATA_DS = DATA_o;
    end
else
    error('ERROR: DATA_o must be 2 or 4 dimensional array')
end