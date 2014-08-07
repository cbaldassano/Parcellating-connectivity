function [D, Pop, Pop_orig, adj_list, borders, names, z_state, FIPS, new_inds] = MergeLowPop(D, Pop, adj_list, contig, borders, names, z_state, FIPS, Pop_thresh)

N = length(Pop);
county_list = cell(N,1);
for i = 1:N
    county_list{i} = i;
end

valid = false(N,1);
valid(contig) = true;

Pop(isnan(Pop)) = 0;
D_nonnorm = diag(Pop)*D;

FIPS_in = FIPS;
FIPS = cell(N,1);
for i = 1:N
    FIPS{i} = {FIPS_in(i,:)};
end
names_in = names;
names = cell(N,1);
for i = 1:N
    names{i} = names_in(i);
end
Pop_orig = cell(N,1);
for i = 1:N
    Pop_orig{i} = {Pop(i)};
end

while (any(Pop(valid) < Pop_thresh))
    [~,min_ind] = min(Pop(valid));
    valid_inds = find(valid);
    min_ind = valid_inds(min_ind);
    
    neighbors = adj_list{min_ind}(ismember(adj_list{min_ind}, valid_inds) & ...
        (z_state(adj_list{min_ind}) == z_state(min_ind))' );
    [~,min_neighbor] = min(Pop(neighbors));
    min_neighbor = neighbors(min_neighbor);
    
    valid(min_neighbor) = false;
    for adj_n = adj_list{min_neighbor}
        if (adj_n ~= min_ind)
            adj_list{adj_n} = unique([adj_list{adj_n} min_ind]);
        end
    end
    adj_list{min_ind} = setdiff([adj_list{min_ind} adj_list{min_neighbor}], min_ind);
    borders{min_ind} = [borders{min_ind};borders{min_neighbor}];
    Pop(min_ind) = Pop(min_ind) + Pop(min_neighbor);
    Pop_orig{min_ind} = [Pop_orig{min_ind}; Pop_orig{min_neighbor}];
    FIPS{min_ind} = [FIPS{min_ind}; FIPS{min_neighbor}];
    names{min_ind} = [names{min_ind}; names{min_neighbor}];
    county_list{min_ind} = [county_list{min_ind} county_list{min_neighbor}];
end


new_inds = find(valid);
new_N = length(new_inds);
new_inds_inv = zeros(N,1);
new_inds_inv(new_inds) = 1:new_N;

D = zeros(new_N);
for i = 1:new_N
    for j = 1:new_N
        if (i ~= j)
            D(i,j) = sum(sum(D_nonnorm(county_list{new_inds(i)},county_list{new_inds(j)})))...
                /Pop(new_inds(i));
        end
    end
end

adj_list_old = adj_list;
adj_list = cell(new_N,1);

for i = 1:new_N
    adj_list{i} = new_inds_inv(adj_list_old{new_inds(i)});
    adj_list{i} = (adj_list{i}(adj_list{i}>0))';
end

Pop = Pop(new_inds);
borders = borders(new_inds);
names = names(new_inds);
FIPS = FIPS(new_inds);
Pop_orig = Pop_orig(new_inds);
z_state = z_state(new_inds);
end