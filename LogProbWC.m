function logp = LogProbWC(D, Z, sizes, alpha, kappa, nu, sigsq)
hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq);

logp = zeros(length(sizes),1);
for i = 1:length(sizes)
    z = cluster(Z, 'maxclust', sizes(i))';
    [sorted_z, sorted_i] = sort(z);
    parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));
    
    % Fake tree c to have correct number of roots
    c = zeros(length(z),1);
    c(1:sizes(i)) = 1:sizes(i);
    
    logp(i) = FullProbabilityddCRP(D, c, parcels, alpha, hyp, CheckSymApprox(D));
    % disp([num2str(sizes(i)) ': ' num2str(logp(i))]);
end