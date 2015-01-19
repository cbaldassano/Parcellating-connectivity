import numpy as np
from scipy import cluster
from scipy import sparse
import random
import ddCRP
import StatsUtil

def InitializeAndRun(Z, D_norm, adj_list, sizes, alpha, kappa, nu, sigsq, pass_limit, gt_z, verbose):

    # Find highest-probability greedy parcellation to initialize
    logp = LogProbWC(D_norm, Z, sizes, alpha, kappa, nu, sigsq)
    max_i = np.argmax(logp)
    z = cluster.hierarchy.fcluster(Z, t=sizes[max_i], criterion = 'maxclust')
  
    # Construct a spanning tree within each cluster as initialization for c
    c = ClusterSpanningTrees(z, adj_list)
    map_z, stats = ddCRP.ddCRP(D_norm, adj_list, c, gt_z, pass_limit, alpha, kappa, nu, sigsq, 1000, verbose)             
    return(map_z, stats)


# Compute probability of each clustering
def LogProbWC(D, Z, sizes, alpha, kappa, nu, sigsq):
    hyp = ddCRP.ComputeCachedLikelihoodTerms(kappa, nu, sigsq)
    
    logp = np.zeros(len(sizes))
    
    for i in range(len(sizes)):
        z = cluster.hierarchy.fcluster(Z, t=sizes[i], criterion = 'maxclust')

        sorted_i = np.argsort(z)
        sorted_z = np.sort(z)
        parcels = np.split(sorted_i,np.flatnonzero(np.diff(sorted_z))+1)
        
        # Fake tree c to have correct number of roots
        c = np.zeros(len(z))
        c[0:sizes[i]] = np.arange(sizes[i])
            
        logp[i] = ddCRP.FullProbabilityddCRP(D, c, parcels, alpha, hyp, StatsUtil.CheckSymApprox(D))

    return logp

    
# Create spanning trees within each cluster, have each element point to its parent
def ClusterSpanningTrees(z, adj_list):

    # We're going to remove edges from adj_list, so make a copy
    adj_list = adj_list.copy()
    
    nvox = len(adj_list)
    # Remove all adjacencies that cross clusters
    for i in range(nvox):
        adj_list[i] = adj_list[i][z[adj_list[i]]==z[i]]
        adj_list[i]  = np.random.permutation(adj_list[i])

    # Construct sparse adjacency matrix
    neighbor_count = [len(neighbors) for neighbors in adj_list]
    node_list = np.zeros(sum(neighbor_count))
    next_edge = 0
    for i in range(nvox):
        if neighbor_count[i] > 0:
            node_list[next_edge:(next_edge+neighbor_count[i])] = i
            next_edge = next_edge + neighbor_count[i]
    G = sparse.csc_matrix((np.ones(len(node_list)),(node_list,np.hstack(adj_list))), shape=(nvox,nvox)) 
    
    # Construct spanning tree in each cluster
    minT = sparse.csgraph.minimum_spanning_tree(G) #creates minimum spanning tree
    c = np.zeros(len(adj_list))
    for clust in np.unique(z):
        clust_vox = np.flatnonzero(z==clust)
        rand_root=clust_vox[random.randint(1,len(clust_vox)-1)]
        _,parents = sparse.csgraph.breadth_first_order(minT,rand_root, directed=False) 
        c[clust_vox] = parents[clust_vox] 
    
    # Roots have parent value of -9999, set them to be their own parent
    roots = np.flatnonzero(c==-9999) 
    c[roots] = roots

    return c
