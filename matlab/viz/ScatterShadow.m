function ScatterShadow(X, radius, Xcol)
scatter3(X(:,1),X(:,2),X(:,3), radius, Xcol, 'filled');
hold on;
N = size(X,1);

xl = get(gca,'XLim');
scatter3(xl(2)*ones(N,1),X(:,2),X(:,3), radius, [0.5 0.5 0.5], 'filled');

yl = get(gca,'YLim');
scatter3(X(:,1),yl(2)*ones(N,1),X(:,3), radius, [0.5 0.5 0.5], 'filled');

zl = get(gca,'ZLim');
scatter3(X(:,1),X(:,2),zl(1)*ones(N,1), radius, [0.5 0.5 0.5], 'filled');
end