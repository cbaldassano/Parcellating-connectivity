function [C_z K varexp] = ClusterSV(C)

sigsqs = logspace(-3, 1, 12);

adj = cell(size(C,1),1);
for i = 1:length(adj)
    adj{i} = [1:(i-1) (i+1):size(C,1)];
end
C_z = zeros(length(sigsqs), size(C,1));
varexp = zeros(length(sigsqs),1);
K = zeros(length(sigsqs),1);
parfor s = 1:length(sigsqs)
    C_z(s,:) = ddCRP(C, adj, [], [], [], 20, 10, 0.0001, 1, sigsqs(s), 100, 0);
    K(s) = length(unique(C_z(s,:)));
    varexp(s) = CalcVarianceExplained(C, C_z(s,:));
end

end