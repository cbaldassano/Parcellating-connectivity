function Z = t_local_linkage(Y,NB,method)

%%
[k, n] = size(Y);
if k~=1
    error('t_local_linkage: ERROR, input must be vector with distnaces.')
end
m = (1+sqrt(1+8*n))/2;
if rem(m,1)
    error('t_local_linkage: ERROR, input must be of length ((2m-1)^2-1)/8, where m is integer.')    
end
   
Z=zeros(m-1,3);

NBM=sparse(size(NB,1));
for ctr=1:size(NB,1)
	nbrs=NB(ctr,:);
	nbrs(nbrs==0)=[];
	NBM(ctr,nbrs)=1;
	NBM(nbrs,ctr)=1;
end
NBM=NBM-diag(diag(NBM));
MASK=squareform(NBM);
%full(squareform(MASK))
%%
% THE FOLLOWING IS "BORROWED" FROM MATLAB ONLY CHANGE IS INTRODUCTION OF MASK

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;


% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
if ~isempty(strmatch(method,['ce';'me';'wa']))
Y = Y .* Y;
end

%%
for s = 1:(n-1)
    fM=find(MASK);
	if strcmp(method,'av')
		p = (m-1):-1:2;
		I = zeros(m*(m-1)/2,1);
		I(cumsum([1 p])) = 1;
		I = cumsum(I);
		J = ones(m*(m-1)/2,1);
		J(cumsum(p)+1) = 2-p;
		J(1)=2;
		J = cumsum(J);
		W = N(R(I)).*N(R(J));
		[v, kk] = min(Y(fM)./W(fM));
	else
		[v, kk] = min(Y(fM));
    end
    if isempty(kk)
        if strcmp(method,'av')
            p = (m-1):-1:2;
            I = zeros(m*(m-1)/2,1);
            I(cumsum([1 p])) = 1;
            I = cumsum(I);
            J = ones(m*(m-1)/2,1);
            J(cumsum(p)+1) = 2-p;
            J(1)=2;
            J = cumsum(J);
            W = N(R(I)).*N(R(J));
            [v, k] = min(Y./W);
        else
            [v, k] = min(Y);
        end
    else 
        k=fM(kk);
    end

    i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
    j = k - (i-1)*(m-i/2)+i; % -> k=j-(i-1)*(m-i/2)+i

    Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A

    % Update Y. In order to vectorize the computation, we need to compute
    % all the indices corresponding to cluster i and j in Y, denoted by I
    % and J.
    I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
    U = [I1 I2 I3];
    I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
    J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

    switch method
        case 'si' % single linkage
        Y(I) = min(Y(I),Y(J));
        case 'co' % complete linkage
        Y(I) = max(Y(I),Y(J));
        case 'av' % average linkage
        Y(I) = Y(I) + Y(J);
        case 'we' % weighted average linkage
        Y(I) = (Y(I) + Y(J))/2;
        case 'ce' % centroid linkage
        K = N(R(i))+N(R(j));
        Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v)./K)./K;
        case 'me' % median linkage
        Y(I) = (Y(I) + Y(J))/2 - v /4;
        case 'wa' % Wards linkage
        Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
        N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
    end
    MASK(I) = MASK(I) | MASK(J); % update mask

    J = [J i*(m-(i+1)/2)-m+j]; % WHAT DOES THIS BIT DO? THIS IS A SINGLE INDEX, PROBABLY THAT ASSOCIATED TO DISTANCE BETWEEN i and j, which obviously does not need updating, but needs deleting!
    Y(J) = [];  % no need for the cluster information about j.
    MASK(J)=[]; % no need for MASK either

    % update m, N, R
    m = m-1;
    N(n+s) = N(R(i)) + N(R(j));
    R(i) = n+s;
    R(j:(n-1))=R((j+1):n);
end

if ~isempty(strmatch(method,['ce';'me';'wa']))
	 Z(:,3) = sqrt(Z(:,3));
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);

