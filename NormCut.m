function z = NormCut(D, adj_list, n_clust, vox_to_clust)
if (nargin < 4)
    vox_to_clust = 1:size(D,1);
end

addpath('normcut');
N = size(D,1);

W = zeros(size(D));
for i = 1:size(D,1)
    for j = adj_list{i}
        W(i,j) = corr2(D(i, (1:N ~= i) & (1:N ~= j)), ...
                       D(j, (1:N ~= i) & (1:N ~= j)));
        if (W(i,j) < 0)
            W(i,j) = 0;
        end
    end
end

nc_discrete = ncutW(W(vox_to_clust,vox_to_clust), n_clust);
z = nc_discrete * (1:n_clust)';
end

