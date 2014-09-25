function D_corr = MapCorrelation(D)
D_corr = zeros(size(D));
valid = true(size(D,1),1);
for i = 1:size(D,1)
    valid(i) = false;
    for j = (i+1):size(D,1)
        valid(j) = false;
        D_corr(i,j) = corr2(D(i,valid),D(j,valid));
        D_corr(j,i) = D_corr(i,j);
        valid(j) = true;
    end
    valid(i) = true;
end