function [X, Y, P, Xpdf, Ypdf, indep, unsmoothP] = CompareConnTypes(conn_vectors)

diff_lim = 5;
n = size(conn_vectors,1);
bin = zeros(n,2,'uint8');
nbins = [100 diff_lim*10];
edges1 = linspace(-0.1,0.5, nbins(1)+1);
X = edges1(1:end-1) + .5*diff(edges1);
edges2 = linspace(0.001, diff_lim, nbins(2)+1);
Y = edges2(1:end-1) + .5*diff(edges2);
Y = cellfun(@num2str,num2cell(round(10*Y)/10),'UniformOutput',false);
Y = ['=0' Y '>5'];
[~,bin(:,1)] = histc(conn_vectors(:,1),edges1);
[~,bin(:,2)] = histc(conn_vectors(:,2),[edges2 Inf]);
bin(:,2) = bin(:,2)+1;
bin = bin(all(bin,2) & (bin(:,1) <= nbins(1)),:);
disp([num2str((n-length(bin))/n*100) '% samples excluded']);
P = accumarray(bin,1,[nbins(1) nbins(2)+2]) ./ n;
unsmoothP = P;

lambda = 1;
P = smooth1D(P,lambda);
P = smooth1D(P',lambda)';

P = P ./sum(P(:));

Xpdf = sum(P,2);
Ypdf = sum(P,1);

indep = Xpdf*Ypdf;

figure;
imagesc(X,1:(nbins(2)+2),(P.*log(P./indep))');
set(gca,'YDir','normal');
set(gca,'YTick',[
set(gca,'YTickLabel',Y);
colormap(PTcolormap(300));
caxis([-1 1]*5*10^(-5));
colorbar;
xlabel('Functional Strength');
ylabel('Anatomical Strength');

% P_x_giv_y = P*diag(1./Ypdf);
% 
% figure;
% KL_giv_y = (P_x_giv_y.*log(P_x_giv_y./repmat(Xpdf,1,size(P,2))))';
% KL_giv_y(isnan(KL_giv_y)) = 0;
% imagesc(X,1:(nbins(2)+2),KL_giv_y);
% set(gca,'YDir','normal');
% set(gca,'YTickLabel',Y);
% colormap(PTcolormap(300));
% caxis([-1 1]* 4*10^-3);
% colorbar;
% xlabel('Functional Strength');
% ylabel('Anatomical Strength');

% figure;
% Xpdf_y0 = P(:,1)/Ypdf(1);
% plot(X, Xpdf_y0.*log(Xpdf_y0./Xpdf));
end

function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
end