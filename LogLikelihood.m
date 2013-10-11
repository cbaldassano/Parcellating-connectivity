function logp = LogLikelihood(stats, hyp)

% stats = [N | mu | sumsq]
% hyp = [mu0 kappa0 nu0 sigsq0 nu0*sigsq0 const_logp_terms]

kappa = hyp(2) + stats(:,1);
nu = hyp(3) + stats(:,1);
%nu_sigsq = hyp(5) + sumSqX + (n*hyp(2)) / (hyp(2)+n) * (hyp(1) - meanX)^2;
%Assume mu0=0 and kappa0 << n
nu_sigsq = hyp(5) + stats(:,3) + hyp(2) * stats(:,2).^2;

logp = hyp(6);
logp = logp + gammaln(nu/2);
logp = logp - 0.5*log(kappa);
logp = logp - (nu/2).*log(nu_sigsq);
logp = logp - (stats(:,1)/2)*log(pi);
logp = sum(logp);
end

