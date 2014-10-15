function cmap = BR2Dcmap(n)

cmap = zeros(n,n,3);
cmap(:,:,3) = 70/100;
for i = 1:n
    for j = 1:n
        cmap(i,j,1) = (238*i+360*j)/(i+j)/360;
        cmap(i,j,2) = max(i,j)/n;
    end
end
cmap = hsv2rgb(cmap);

end