function var_exp = CalcVarianceExplained(D, z, base_var)

if (iscell(z))
    bins = z;
else
    z = z(:)';
    [sorted_z, sorted_i] = sort(z);
    bins = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));
end

D = cast(D,'double');
if (nargin < 3)
    off_diags = true(size(D));
    for i = 1:size(D,1)
        off_diags(i,i) = false;
    end
    off_diags = off_diags(:);
    base_var = sum((mean(D(off_diags)) - D(off_diags)).^2);
end
sum_var = 0;
for i=1:length(bins)
    %fprintf(' %d', i);
    for j=i:length(bins)
        x = D(bins{i},bins{j});
        if (i==j)
            x = x(~tril(ones(size(x))));
        else
            x = x(:);
        end
        sum_var = sum_var + sum((x - mean(x)).^2);
    end
end

var_exp = 1 - 2 *sum_var / base_var;
end

