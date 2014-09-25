function D = LoadTractography(nearest_vertex, probtrackx_path)

D = sparse(59412,59412);
for f = 1:20
    disp(['Loading streamlines ' num2str(f) '/20...']);
    matrix3 = importdata([probtrackx_path num2str(f) '/fdt_matrix3.dot']);
    matrix3(:,1:2) = nearest_vertex(matrix3(:,1:2));
    matrix3 = matrix3(all(matrix3(:,1:2)>0,2),:);
    
    D = D + sparse(matrix3(:,1),matrix3(:,2),matrix3(:,3),59412,59412);
end
D = single(full(D));
for i = 1:59412
    D(i,i) = 0;
end
D = D + D';
end