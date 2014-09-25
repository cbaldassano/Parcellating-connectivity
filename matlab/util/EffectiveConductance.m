function G = EffectiveConductance(M)

if (M(1,1) == 0)
    % Assume this is a weight matrix
    Lpinv = pinv(diag(sum(M)) - M);
else
    % Assume this the pseudo-inv of the laplacian
    Lpinv = M;
end

N = size(Lpinv,1);
Lxx = repmat(diag(Lpinv),1,N);
G = 1./(Lxx + Lxx' - 2*Lpinv);

end