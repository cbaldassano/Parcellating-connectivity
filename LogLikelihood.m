function logp = LogLikelihood(x, hyp)

if (isempty(x))
    logp = 0;
    return;
end

mu0 = 0;
kappa0 = hyp(1);
nu0 = hyp(2);
sigsq0 = hyp(3);

n = length(x);
meanX = sum(x)/n;
sumSqX = sum((x-meanX).^2);

kappa = kappa0 + n;
nu = nu0 + n;
sigsq = (1/nu) * (nu0*sigsq0 + sumSqX + ...
                 (n*kappa0) / (kappa0 + n) * (mu0 - meanX)^2);

logp = gammaln(nu/2)-gammaln(nu0/2) + ...
    (1/2)*log(kappa0) - (1/2)*log(kappa) +...
    (nu0/2)*log(nu0*sigsq0) - (nu/2)*log(nu*sigsq) - (n/2)*log(pi);

end

