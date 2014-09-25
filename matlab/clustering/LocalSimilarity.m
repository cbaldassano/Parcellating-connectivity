function z = LocalSimilarity(D, adj_list, n_clust, vox_to_clust)
if (nargin < 4)
    vox_to_clust = 1:size(D,1);
end

N = size(D,1);
W = inf(size(D));
if (CheckSymApprox(D))
    for i = 1:size(D,1)
        for j = adj_list{i}
            W(i,j) = norm(D(i, (1:N ~= i) & (1:N ~= j)) - ...
                          D(j, (1:N ~= i) & (1:N ~= j)), 2);
        end
    end
else
    for i = 1:size(D,1)
        for j = adj_list{i}
            W(i,j) = norm([D(i, (1:N ~= i) & (1:N ~= j)) D((1:N ~= i) & (1:N ~= j), i)'] - ...
                          [D(j, (1:N ~= i) & (1:N ~= j)) D((1:N ~= i) & (1:N ~= j), j)'], 2);
        end
    end
end


z = ThresholdSimilarity(W(vox_to_clust,vox_to_clust), n_clust);
end

function z = ThresholdSimilarity(W, n_clust)

thresh = mean(W(isfinite(W)));
thresh_step = 0.1;
last_increase = 0;
last_decrease = 0;
while(1)
    G = sparse(W <= thresh);
    [~, z] = graphconncomp(G);
    curr_n_clust = length(unique(z));
    if (curr_n_clust == n_clust)
        break;
    elseif (curr_n_clust < n_clust)
        if (last_increase)
            thresh_step = thresh_step^2;
        end
        thresh = thresh * (1-thresh_step);
        last_decrease = 1; last_increase = 0;
    else
        if (last_decrease)
            thresh_step = thresh_step^2;
        end
        thresh = thresh * (1+thresh_step);
        last_decrease = 0; last_increase = 1;
    end
end

end