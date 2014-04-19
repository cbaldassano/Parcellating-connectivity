function [sq_diff_map] = ConnDiffMap(conn_vectors)

reg_subset = randperm(size(conn_vectors,1),10^6);
b = RegressSimple(conn_vectors(reg_subset,1),conn_vectors(reg_subset,2));
mean_func = mean(conn_vectors(reg_subset,1));

N = 59412;
sq_diff_map = zeros(N,1);

for v = 1:N
    if (mod(v,1000)==1)
        disp([num2str(v/N*100) '%...']);
    end
    inds = [((v^2-v)/2-max(v-2,0)):(v^2-v)/2 ((v^2-v)/2+1+cumsum((v-1):(N-2)))];
    inds = inds(inds>0);
    sq_diff_map(v) = sum((conn_vectors(inds,1) - mean_func).^2) - ...
                     sum((conn_vectors(inds,1) - b(1) - b(2)*conn_vectors(inds,2)).^2);
end