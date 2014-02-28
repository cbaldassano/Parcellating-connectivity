function grad_map = GradientMap(Y, adj_list)

if (size(Y,1) > 1)
    Qy = repmat(dot(Y,Y,2),1,size(Y,1));
    Y = Qy+Qy'-2*(Y*Y');
    Y(1:(size(Y,1)+1):size(Y,1)^2) = 0; % Remove numerical errors on diagonal
    Y = squareform(Y);
end

n = length(adj_list);
grad_map = zeros(n,1);
valid_dims = true(n,1);
for i = 1:n
    if (mod(i,1000)==0)
        disp(['V ' num2str(i)]);
    end
    local_grad = zeros(length(adj_list{i}),1);
    valid_dims(i) = 0;
    for j = 1:length(adj_list{i})
        valid_dims(adj_list{i}(j)) = 0;
        local_grad(j) = norm(Y(PdistInds(i, n, valid_dims)) - Y(PdistInds(adj_list{i}(j), n, valid_dims)), 2)^2;
        valid_dims(adj_list{i}(j)) = 1;
    end
    valid_dims(i) = 1;
    grad_map(i) = mean(local_grad);
end

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