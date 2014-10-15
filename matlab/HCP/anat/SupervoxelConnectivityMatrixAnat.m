function C = SupervoxelConnectivityMatrixAnat(D, z, alpha)

[sorted_z, sorted_i] = sort(z);
parcels = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));

C = zeros(length(parcels),length(parcels));

%disp('Computing matrix...');
for i = 1:length(parcels)
    %fprintf('%d ', i);
    for j = (i+1):length(parcels)
         C(i,j) = log(mean(mean(exp(D(parcels{i},parcels{j})))));
%         if (nargin < 3 || alpha == 0)
%             C(i,j) = mean(mean(D(parcels{i},parcels{j})));
%         elseif (isinf(alpha))
%             C(i,j) = max(max(D(parcels{i},parcels{j})));
%         else
%             vals = D(parcels{i},parcels{j});
%             vals = vals(:);
%             weights = exp(alpha*(vals-max(vals)));
%             C(i,j) = sum(vals.*weights)/sum(weights);
%         end
            
    end
end

C = C + C';


end

