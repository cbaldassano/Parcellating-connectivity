function [map_z, stats] = InitializeAndRunddCRP(Z, D_norm, adj_list, sizes, alpha, kappa, nu, sigsq, pass_limit, gt_z, verbose)
% Standard alpha = 10, kappa = 0.0001, nu = 1

logp = LogProbWC(D_norm, Z, sizes, alpha, kappa, nu, sigsq);
[~,max_i] = max(logp);
z = cluster(Z, 'maxclust', sizes(max_i));
c = ClusterSpanningTrees(z, adj_list);
[map_z,stats] = ddCRP(D_norm, adj_list, c, [], gt_z, ...
                  pass_limit, alpha, kappa, nu, sigsq, ...
                  1000, verbose);
end