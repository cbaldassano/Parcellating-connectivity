function logp = FullProbabilityddCRP(D, c, parcels, alpha, hyp)

logp = log(alpha) * sum(c == 1:length(c));
for c1 = 1:length(parcels)
    for c2 = c1:length(parcels)
        samples = D(parcels{c1},parcels{c2});
        if (c1 == c2)
            samples = samples(logical(triu(ones(size(samples)),1)));
        else
            samples = samples(:);
        end
        logp = logp + LogLikelihood(samples, hyp);
    end
end

end

