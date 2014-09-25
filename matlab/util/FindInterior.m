function FindInterior(input_nii, xrange, output_nii)

inp = load_nii(input_nii);
yrange = 1:size(inp.img,2);
zrange = 1:size(inp.img,3);
xvol = inp.img;
yvol = inp.img;
zvol = inp.img;

for y = 1:size(inp.img,2)
    for z= 1:size(inp.img,3)
        xvol(xrange,y,z) = mod([0;cumsum(max(diff(inp.img(xrange,y,z)),0))],2);
    end
end

for x = 1:size(inp.img,1)
    for z= 1:size(inp.img,3)
        yvol(x,yrange,z) = mod([0 cumsum(max(diff(squeeze(inp.img(x,yrange,z))),0))],2);
    end
end

for x = 1:size(inp.img,1)
    for y= 1:size(inp.img,2)
        zvol(x,y,zrange) = mod([0;cumsum(max(diff(squeeze(inp.img(x,y,zrange))),0))],2);
    end
end

inp.img = inp.img + ((xvol + yvol + zvol) == 3);
save_nii(inp, output_nii);

end

