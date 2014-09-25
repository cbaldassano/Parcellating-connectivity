function [b R2 r] = RegressSimple(y,X)
X_c = [ones(size(X,1),1) X];
[Q,R] = qr(X_c,0);
b = R \ (Q' * y);

if (nargout > 1)
    yhat = X_c*b;
    r = y - yhat;
    R2 = zeros(size(y,2),1);
    for i = 1:size(R2,1)
        R2(i) = 1 - sum(r(:,i).^2)/sum((y(:,i)-mean(y(:,i))).^2);
    end
end
end