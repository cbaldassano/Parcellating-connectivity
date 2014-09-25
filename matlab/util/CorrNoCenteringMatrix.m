function R = CorrNoCenteringMatrix(V)

norms = sqrt(sum(V.^2, 2));
V = diag(1./norms)*V;
R = V*V';

end

