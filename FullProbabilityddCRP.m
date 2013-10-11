function logp = FullProbabilityddCRP(D, c, parcels, alpha, hyp)

stats = zeros(length(parcels)*(length(parcels)+1),3);
j = 1;
for c1 = 1:length(parcels)
    for c2 = c1:length(parcels)
        samples = D(parcels{c1},parcels{c2});
        if (c1 == c2)
            samples = samples(logical(triu(ones(size(samples)),1)));
            if (isempty(samples))
                continue;
            end
        else
            samples = samples(:);
        end
        stats(j,1) = length(samples);
        stats(j,2) = sum(samples)/stats(j,1);
        stats(j,3) = sum((samples-stats(j,2)).^2);
        j = j+1;
    end
end
logp = log(alpha) * sum(c' == 1:length(c)) + LogLikelihood(stats, hyp);

end

