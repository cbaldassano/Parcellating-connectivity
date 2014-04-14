function var_exp = VarExpTree(D, Z, sizes, base_var)

    D = cast(D, 'double');
    if (nargin < 4)
        off_diags = true(size(D));
        for i = 1:size(D,1)
            off_diags(i,i) = false;
        end
        off_diags = off_diags(:);
        base_var = sum((mean(D(off_diags)) - D(off_diags)).^2);
    end

    var_exp = zeros(length(sizes),1);
    
    for i = 1:length(sizes)
        z = cluster(Z, 'maxclust', sizes(i));
        var_exp(i) = CalcVarianceExplained(D, z', base_var);
        disp([num2str(sizes(i)) ': ' num2str(var_exp(i))]);
    end
end
