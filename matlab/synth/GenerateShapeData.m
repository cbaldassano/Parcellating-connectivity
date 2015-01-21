function [D adj_list z coords] = GenerateShapeData(sig)

data = importdata('synth/spiral.txt');
coords = data(:,1:2);
z = data(:,3);

N = length(z);
K = length(unique(z));

A = randn(K,K);

D = zeros(N, N);
for v1 = 1:N
    for v2 = 1:N
        if (v1 ~= v2)
            D(v1,v2) = sig*randn(1) + A(z(v1),z(v2));
        end
    end
end

adj_list = cell(303,1);
for i = 1:303
    adj_list{i} = [];
    if (mod(i,101) ~= 1)
        adj_list{i} = [adj_list{i} i-1];
    end
    if (mod(i,101) ~= 0)
        adj_list{i} = [adj_list{i} i+1];
    end
end

for i = 33:101
    adj_list{i} = [adj_list{i} i+170];
    adj_list{i+170} = [adj_list{i+170} i];
end
for i = 133:202
    adj_list{i} = [adj_list{i} i-132];
    adj_list{i-132} = [adj_list{i-132} i];
end
for i = 236:303
    adj_list{i} = [adj_list{i} i-134];
    adj_list{i-134} = [adj_list{i-134} i];
end
    

end

