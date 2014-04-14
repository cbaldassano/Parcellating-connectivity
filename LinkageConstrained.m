function Z = LinkageConstrained(X, adj_list)
% Based on MATLAB "LINKAGEOLD" function (in linkage.m)

% load('../data/Q3/group/smoothed1_sqEucDist.mat');
Qx = repmat(dot(X,X,2),1,size(X,1));
Y = Qx+Qx'-2*(X*X');
Y(1:(size(Y,1)+1):size(Y,1)^2) = 0; % Remove numerical errors on diagonal
Y = squareform(Y);

A = false(length(adj_list));
for i = 1:length(adj_list)
    A(i,adj_list{i}) = true;
end
connected = squareform(A);

N = length(adj_list);
valid_clusts = true(N,1);
col_limits = cumsum((N-1):-1:1);

m = ceil(sqrt(2*size(Y,2))); % (1+sqrt(1+8*n))/2, but works for large n
if isa(Y,'single')
   Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
   Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  C denotes how many points are contained in each cluster.
C = zeros(1,2*m-1);
C(1:m) = 1;
R = 1:m;

all_inds = 1:size(Y,2);
conn_inds = all_inds(connected);

for s = 1:(m-1)
   %fprintf('%d%% clustered...\n', round(100*s/(m-1)));
   
   if (isempty(conn_inds))
       % The graph was disconnected (e.g. two hemispheres)
       % Just add all connections to finish up cluster tree
       connected = zeros(length(connected),1);
       conn_inds = [];
       valid_clust_inds = find(valid_clusts);
       for i = valid_clust_inds'
           U = valid_clusts;
           U(i) = 0;
           new_conns = PdistInds(i, N, U);
           connected(new_conns) = 1;
           conn_inds = [conn_inds new_conns];
       end
       conn_inds = unique(conn_inds);
   end
   
   [v, k] = min(Y(conn_inds));
   k = conn_inds(k);
   
   j = find(k <= col_limits, 1, 'first');
   i = N - (col_limits(j) - k);

   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A

   % Update Y. In order to vectorize the computation, we need to compute
   % all the indices corresponding to cluster i and j in Y, denoted by I
   % and J.
   U = valid_clusts;
   U([i j]) = 0;
   I = PdistInds(i, N, U);
   J = PdistInds(j, N, U);

   Y(I) = ((C(R(U))+C(R(i))).*Y(I) + (C(R(U))+C(R(j))).*Y(J) - ...
    C(R(U))*v)./(C(R(i))+C(R(j))+C(R(U)));

   % update connected
   new_conns = connected(J) & ~connected(I);
   connected(I) = connected(I) | new_conns;
   conn_inds = sort([conn_inds I(new_conns)]);
   
   U(i) = 1;
   J = PdistInds(j, N, U);
   conn_inds(ismembc(conn_inds,J)) = [];
   connected(J) = false(1, length(J));
   
   valid_clusts(j) = 0;
   
   % update m, N, R
   C(m+s) = C(R(i)) + C(R(j));
   R(i) = m+s;
end

Z(:,3) = sqrt(Z(:,3));

end


function I = PdistInds(row, N, valid_flags)

if (row > 1)
    inds1 =[(row-1) (row-1)+cumsum((N-2):-1:(N-row+1))];
    I = [inds1 0 (inds1(end)+N-row+1):(inds1(end)+2*N-2*row)];
else
    I = 0:(N-1);
end

I = I(valid_flags);
    
end
