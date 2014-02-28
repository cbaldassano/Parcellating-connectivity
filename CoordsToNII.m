function CoordsToNII(coords, vals, max_dist, ref_file, out_file)

if (length(vals) == 1)
    vals = vals*ones(size(coords,1),1);
end

ref = load_nii(ref_file);
ref_dim = ref.hdr.dime.dim(2:4);
Smat = [ref.hdr.hist.srow_x; ref.hdr.hist.srow_y; ref.hdr.hist.srow_z];

if (Smat(1,1) < 0)
    Smat(1,1) = abs(Smat(1,1));
    Smat(1,4) = Smat(1,4) - (ref_dim(1)-1)*Smat(1,1);
    ref.hdr.hist.qoffset_x = Smat(1,4);
end
ref.img = zeros(ref_dim(1), ref_dim(2), ref_dim(3), 'single');
vox_scales = diag(Smat)';
radius = floor(max(max_dist./vox_scales));

for i = 1:size(coords,1)
    center = MMToInd(Smat, coords(i,:));
    for x = -radius:radius
        for y = -radius:radius
            for z = -radius:radius
                if (norm([x y z].*vox_scales) <= max_dist)
                    ref.img(center(1)+x, center(2)+y, center(3)+z) = vals(i);
                end
            end
        end
    end
end

save_nii(ref, out_file);

end

function ind = MMToInd(Smat, mm)
    ind = round(Smat(1:3,1:3)\(mm(:) - Smat(:,4)));
end