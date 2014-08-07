function [z_size] = ShowConnMat(D_nonnorm, Pop, z)
% t = tabulate(z);
% t = t(:,2);
% [~,sizeorder] = sort(t,'descend');
% z_size = zeros(length(z),1);
% for i = 1:max(z)
%     z_size(z==sizeorder(i)) = i;
% end

clustPop = zeros(max(z),1);
for i = 1:max(z)
    clustPop(i) = sum(Pop(z==i));
end
[~,Poporder] = sort(clustPop,'descend');
z_size = zeros(length(z),1);
for i = 1:max(z)
    z_size(z==Poporder(i)) = i;
end

[~,sorted] = sort(z_size,'ascend');

total_movers = sum(sum(D_nonnorm));

figure;
imagesc(D_nonnorm(sorted(1:1051),sorted(1:1051))*sum(Pop)^2/total_movers.*repmat(1./Pop(sorted(1:1051)),1,1051).*repmat(1./Pop(sorted(1:1051))',1051,1))
caxis([0 10]); colormap(PTcolormap(200,[-1 9]))

C = zeros(max(z_size));
%C_nonnorm = zeros(max(z_size));
%pmat = NaN(max(z_size));
for i = 1:max(z_size)
    for j = 1:max(z_size)
        C(i,j) = sum(sum(D_nonnorm(z_size==i, z_size==j)))/(sum(Pop(z_size==i))*sum(Pop(z_size==j)));
        %C_nonnorm(i,j) = sum(sum(D_nonnorm(z_size==i, z_size==j)));
        %pmat(i,j) = 1-binocdf(C_nonnorm(i,j),total_movers,sum(Pop(z_size==i))*sum(Pop(z_size==j))/(sum(Pop)^2));
    end
end
C_norm = C(1:10,1:10)*sum(Pop)^2/total_movers;
figure; imagesc(C_norm);
caxis([0 10]); colormap(PTcolormap(200,[-1 9]));
[r,c] = find(C_norm > 1);
for i = 1:length(r)
    disp([num2str(r(i)) ' to ' num2str(c(i)) ': ' num2str(C_norm(r(i),c(i)))]);
end
%imagesc(C(z_size(sorted),z_size(sorted)));


end