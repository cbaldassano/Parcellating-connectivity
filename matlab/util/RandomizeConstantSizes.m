function z = RandomizeConstantSizes(z, adj_list, restarts)

s_target = ParcelSizes(z);
s = s_target;
z_orig = z;
for r = 1:restarts
    if (r==restarts || mod(r,5)==0)
        disp(['Restart ' num2str(r) ', NMI=' num2str(CalcNMI(z_orig,z))]);
    end
    for i = 1:max(z)
        [z s] = AddNeighbor(z, adj_list, s, i);
    end
    fprintf('  Missized parcels:     ');
    if (r < restarts)
        missized_target = 1;
    else
        missized_target = 0;
    end
    while (sum(s < s_target) > missized_target)
        small = find(s < s_target);
        fprintf('\b\b\b\b%4d', length(small));
        [z s] = AddNeighbor(z, adj_list, s, small(randi(length(small))));
    end 
    fprintf('\n');
end

end

function [z s] = AddNeighbor(z, adj_list, s, z_to_add)

neighbors = setdiff(horzcat(adj_list{z==z_to_add}), find(z==z_to_add));
neighbors = setdiff(neighbors, find(ismember(z,find(s==1))));
if (~isempty(neighbors))
    z(neighbors(randi(length(neighbors)))) = z_to_add;
    s = ParcelSizes(z);
end

end

function s = ParcelSizes(z)
    s = tabulate(z);
    s = s(:,2);
end