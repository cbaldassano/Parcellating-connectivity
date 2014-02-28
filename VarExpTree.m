function var_exp = VarExpTree(D, Z, sizes, base_var)

    if (nargin < 4)
        off_diags = ~logical(eye(size(D))); off_diags = off_diags(:);
        base_var = sum((mean(D(off_diags)) - D(off_diags)).^2);
    end

    n = size(D,1);
    var_exp = zeros(length(sizes),1);
    
    parcels = [num2cell(1:n) cell(1,n-1)];
    for i = 1:size(Z,1)
        j = find(sizes == (n - i + 1));
        if (~isempty(j))
            disp(num2str(sizes(j)));
            parcels_stripped = parcels(~cellfun('isempty', parcels));
            var_exp(j) = CalcVarianceExplained(D, parcels_stripped, base_var);
        end
        parcels{n+i} = [parcels{Z(i,1)} parcels{Z(i,2)}];
        parcels{Z(i,1)} = [];
        parcels{Z(i,2)} = [];
    end
end