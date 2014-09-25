function WriteCoordsBrik(templateFile, coords, voxHalfSide, vals, saveFile)

[err, V, Info, ErrMessage] = BrikLoad(templateFile);
origin = Info.ORIGIN;
step_size = Info.DELTA;

coords = round(bsxfun(@minus, coords, origin)*diag(1./step_size));
V = zeros(size(V));
for i = 1:size(coords,1)
    for x = -voxHalfSide:voxHalfSide
        for y = -voxHalfSide:voxHalfSide
            for z = -voxHalfSide:voxHalfSide
                V(coords(i,1)+x,coords(i,2)+y,coords(i,3)+z) = vals(i);
            end
        end
    end
end

Opt.Prefix = saveFile;
WriteBrik(V,Info,Opt);
end