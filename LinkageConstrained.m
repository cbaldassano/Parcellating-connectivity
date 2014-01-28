function Z = LinkageConstrained(X, adj_list)
% Based on MATLAB "LINKAGEOLD" function (in linkage.m)

Y = pdist(X,'euclidean');
A = inf(size(X,1));
for i = 1:length(adj_list)
    A(i,adj_list{i}) = 1;
end
A(1:(size(A,1)+1):(size(A,1)^2)) = 0;
connected = squareform(A);

n = size(Y,2);
m = ceil(sqrt(2*n)); % (1+sqrt(1+8*n))/2, but works for large n
if isa(Y,'single')
   Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
   Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;

% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
Y = Y .* Y;

for s = 1:(n-1)
   [v, k] = min(Y.*connected);

   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
   j = k - (i-1)*(m-i/2)+i;

   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A

   % Update Y. In order to vectorize the computation, we need to compute
   % all the indices corresponding to cluster i and j in Y, denoted by I
   % and J.
   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
   U = [I1 I2 I3];
   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

   Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
    N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));

   % update connected
   new_connI = Inf(1,length(I));
   new_connI(isfinite(connected(I)) | isfinite(connected(J))) = 1;
   connected(I) = new_connI;
   
   % no need for the cluster information about j
   J = [J i*(m-(i+1)/2)-m+j];
   Y(J) = [];
   connected(J) = [];
   
   % update m, N, R
   m = m-1;
   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(n-1))=R((j+1):n);
end

Z(:,3) = sqrt(Z(:,3));

Z(:,[1 2])=sort(Z(:,[1 2]),2);
