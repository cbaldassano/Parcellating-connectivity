function D = LoadTractography(nearest_vertex)

D = sparse(59412,59412);
for f = 1:20
    f
    matrix3 = importdata(['/mnt/chekov_scratch/HCP/101915/T1w/probtrackx/' num2str(f) '/fdt_matrix3.dot']);
    matrix3(:,1:2) = nearest_vertex(matrix3(:,1:2));
    matrix3 = matrix3(all(matrix3(:,1:2)>0,2),:);
    
    D = D + sparse(matrix3(:,1),matrix3(:,2),matrix3(:,3),59412,59412);
end
D = single(full(D));
D(logical(eye(size(D)))) = 0;
D = D + D';
end