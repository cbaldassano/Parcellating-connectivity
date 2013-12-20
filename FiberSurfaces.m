function FiberSurfaces(diff_dump_file, pos_surface_file, neg_surface_file)

pts = importdata(diff_dump_file);
pts(:,1:3) = pts(:,1:3).*repmat([-1 -1 1], size(pts,1), 1);
min_pts = round(min(pts(:,1:3))*10)/10;
max_pts = round(max(pts(:,1:3))*10)/10;
dims = round((max_pts - min_pts)/0.7) + 1;
V = zeros(dims(1),dims(2),dims(3));

for i = 1:size(pts,1)
    ind = round((pts(i,1:3) - min_pts)/0.7)+1;
    V(ind(1),ind(2),ind(3)) = pts(i,4);
end

[x,y,z] = meshgrid(min_pts(1):0.7:max_pts(1),min_pts(2):0.7:max_pts(2),min_pts(3):0.7:max_pts(3));
pos_fv = isosurface(x,y,z,permute(V,[2 1 3]),0.1);
neg_fv = isosurface(x,y,z,permute(V,[2 1 3]),-0.1);

stlwrite(pos_surface_file,pos_fv,'mode','binary','invertnormals',true);
stlwrite(neg_surface_file,neg_fv,'mode','binary');

end

