function orig = ConvertToOrigHem(map, orig_ind)

orig = cell(2,1);
orig{1} = zeros(length(orig_ind{1}),1);
orig{2} = zeros(length(orig_ind{2}),1);
orig{1}(orig_ind{1}) = map(1:length(orig_ind{1}));
orig{2}(orig_ind{2}) = map((length(orig_ind{1})+1):end);

end