function t = ParcelTscores(z, data)

[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

t = zeros(length(parcels),1);
for i = 1:length(parcels)
    samp = data(parcels{i});
    samp = samp(isfinite(samp));
    if (length(samp) > 1)
        t(i) = mean(samp)/(std(samp)/sqrt(length(samp)));
    else
        t(i) = 0;
    end
end

end