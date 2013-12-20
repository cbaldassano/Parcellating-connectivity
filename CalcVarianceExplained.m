function var_exp = CalcVarianceExplained(D, z)

[sorted_z, sorted_i] = sort(z);
bins = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

est_D = zeros(size(D));
for i=1:length(bins)
    disp(num2str(i));
    for j=1:length(bins)
        x = D(bins{i},bins{j});
        if (i==j)
            x = x(~tril(ones(size(x))));
        else
            x = x(:);
        end
        est_D(bins{i},bins{j}) = mean(x);
    end
end

off_diags = ~eye(size(D)); off_diags = off_diags(:);
var_exp = 1 - sum((est_D(off_diags) - D(off_diags)).^2) / ...
              sum((mean(D(off_diags)) - D(off_diags)).^2);
end

