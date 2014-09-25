function pts = ImportSavedPaths(saved_path_file)

pts = cell(0,1);
fid = fopen(saved_path_file, 'r');

tract_i = 1;
while(~feof(fid))
    tract_size = fscanf(fid, '# %d\n', 1);
    points = reshape(fscanf(fid, '%f %f %f\n'), 3, tract_size)';
    if (mod(tract_i,3)==0)
        seed = find((points(:,1) == points(1,1)) & ...
                    (points(:,2) == points(1,2)) & ...
                    (points(:,3) == points(1,3)));
        points = [points((seed(2)-1):-1:2,:);points(seed(2):end,:)];
        pts = [pts; cell(1,1)];
        pts{end} = points;
    end
    tract_i = tract_i + 1;
end

fclose(fid);


end

