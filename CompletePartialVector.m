function D_complete = CompletePartialVector(D, N, partial_to_complete)

D_complete = zeros((N^2-N)/2,1);
N_partial = length(partial_to_complete);

for v = 1:N_partial
    partial_inds =[((v^2-v)/2-max(v-2,0)):(v^2-v)/2 ((v^2-v)/2+1+cumsum((v-1):(N_partial-2)))];
    partial_inds = partial_inds(partial_inds > 0);
    partial_vec = D(partial_inds);
    complete_vec = zeros(N-1,1);
    complete_vec([partial_to_complete(1:(v-1)) partial_to_complete((v+1):N_partial)-1]) = partial_vec;
    v2 = partial_to_complete(v);
    complete_inds = [((v2^2-v2)/2-max(v2-2,0)):(v2^2-v2)/2 ((v2^2-v2)/2+1+cumsum((v2-1):(N-2)))];
    complete_inds = complete_inds(complete_inds > 0);
    D_complete(complete_inds) = complete_vec;
end
end