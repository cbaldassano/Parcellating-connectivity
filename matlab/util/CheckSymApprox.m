function sym = CheckSymApprox(D)

sym_sub = [randi(size(D,1), 1000,1) randi(size(D,1), 1000,1)];
sym = all(D(sub2ind(size(D), sym_sub(:,1), sym_sub(:,2)))==...
          D(sub2ind(size(D), sym_sub(:,2), sym_sub(:,1))));
end