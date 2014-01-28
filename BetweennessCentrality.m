function BC = BetweennessCentrality(C, z)

[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

C_sim = C;
C_sim(logical(eye(size(C)))) = 1;
C_sim(C_sim < 0) = 0;
C_inv_e = exp(-1*C_sim)-exp(-1);

BC = betweenness_wei(C_inv_e,cellfun(@length,parcels));
end