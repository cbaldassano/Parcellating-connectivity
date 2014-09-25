function D_sm = SmoothConnMat(D, adj_list)

D_sm = D;
for v1 = 1:size(D,1)
    if (mod(v1,100)==0)
        fprintf('%d ', v1);
    end
    for v2 = (v1+1):size(D,1)
        D_sm(v1,v2) = mean([D(v1,v2);D(v1, adj_list{v2})';D(adj_list{v1}, v2)]);
    end
end
D_sm(tril(true(size(D)))) = 0;
D_sm = D_sm + D_sm';
end