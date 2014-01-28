function WriteMatrix32(M)
    M = cast(M,'single');
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            fprintf('%f, ', M(i,j));
        end
        fprintf('\n');
    end
end