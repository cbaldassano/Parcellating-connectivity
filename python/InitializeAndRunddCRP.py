import numpy as np
import scipy as sp
import random
import math as mt

def InitializeAndRunddCRP(Z, D_norm, adj_list, sizes, alpha, kappa, nu, sigsq, pass_limit, burn_in_passes, gt_z, verbose):
    # Standard alpha = 10, kappa = 0.0001, nu = 1
    logp = LogProbWC(D_norm, Z, sizes, alpha, kappa, nu, sigsq)
    max_i = np.argmax(logp) #same as [~,max_i] = max(logp)...ie the index for max value of vector logp 
    z= sp.cluster.hierarchy.fcluster(Z, t=sizes[max_i], criterion = 'maxclust')  #should be same as z = cluster(Z, 'maxclust', sizes(max_i))    
  
    c = ClusterSpanningTrees(z, adj_list)
    map_z, stats, pair_prob = ddCRP(D_norm, adj_list, c, [], gt_z, pass_limit, alpha, kappa, nu, sigsq, burn_in_passes, 1000, verbose)             
    return(map_z, stats, pair_prob)


def LogProbWC(D_norm, Z, sizes, alpha, kappa, nu, sigsq):
    hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq)
    
    logp = np.zeros(len(sizes))
    
    for i in range(len(sizes)):
        z= sp.cluster.hierarchy.fcluster(Z, t=sizes[i], criterion = 'maxclust') #should be same as z = cluster(Z, 'maxclust', sizes(i))
        sorted_z = np.sort(z, axis=0) #to sort columns like matlab, use axis = 0
        sorted_i = np.argsort(z, axis=0) #the indices when sorted
        
        #following lines creates [0 sorted_z (max(z)+1)]
        diff_input = np.append(sorted_z,np.amax(z)+1)
        diff_input = np.insert(diff_input,0,0)
        
        find_diff = np.nonzero(np.diff(diff_input)) #find(diff([0 sorted_z (max(z)+1)]))
        diff1 = np.diff(find_diff) #diff(find(diff([0 sorted_z (max(z)+1)])))
        
        #the following does equivalent of mat2cell, except creates list of list
        #try the following:
        #sorted_i = np.array([1, 2, 3, 4, 0, 0, 0, 0, 1, 2])
        #diff1 = np.array([2, 2, 2, 4])
        #print parcels and see that sort_i is clustered into lists of size dictated by diff1, e.g. parcels =[[1,2],[3,4],[0,0],[0,0,1,2]]
        parcels = []
        count = 0         
        for j in range(len(diff1)):
            parcels.append(sorted_i[count:count+diff1[j]])
            count = count+diff1[j]
            #print(parcels[j])
        
        c = np.zeros(len(z))
        c[0:sizes[i]] = np.arange(sizes[i]) #c(1:sizes(i)) = 1:sizes(i)
            
        logp[i] = FullProbabilityddCRP(D, c, parcels, alpha, hyp, CheckSymApprox(D))

    return logp

    
    
def ClusterSpanningTrees(z, adj_list):    
    nvox = len(adj_list)
    for i in range(nvox):
        #Unclear what is meant here...if z(adj_list[i]) == z[i], then select value from adj_list[i]
        adj_list[i] = adj_list[i]*(z(adj_list[i]) == z[i]) 
        adj_list[i]  = adj_list[i] * (np.random.permutation(range(len(adj_list[i])))) #unclear what is being done here..
    neighbor_count = [len(neighbors) for neighbors in adj_list] #same as cellfun(@length, adj_list); implemented in python as a list comprehension that counts members in each list (within the list)
    node_list = np.zeros(sum(neighbor_count)) #neighbor_count is list type
    next_edge = 0 #same as next_edge = 1?
    
    for i in range(nvox):
        if neighbor_count[i] > 0:
            node_list[next_edge:(next_edge+neighbor_count[i]-1)] = i
            next_edge = next_edge + neighbor_count[i]
    #csc_matrix( (data,(row,col)), shape=(n,n)) 
    #same as G = sparse(node_list, [adj_list{:}]', 1, nvox, nvox). Can't use 1, must use ones(nvox)
    G = sp.sparse.csc_matrix((np.ones(nvox),(np.arange(nvox),np.array(adj_list))), shape=(nvox,nvox)) 
    
    c = np.zeros(len(adj_list))
    for clust in np.unique(z):
        clust_vox = np.nonzero(z==clust) #clust_vox = find(z==clust).  Finds all indices in z such that that all have same clust value
        clust_vox=clust_vox.tolist() #convert clust_vox to list because later parents[clust_vox] works only if clust_vox is type list (or iterable)
        
        #following is [~,parents] = graphminspantree(G,clust_vox(randi(length(clust_vox),1)))
        rand_root=clust_vox[random.randint(1,len(clust_vox))]
        minT = sp.sparse.csgraph.minimum_spanning_tree(G) #creates minimum spanning tree
        #parents is a vector with root node's parent assigned -9999.  To create paretns, we use breath-first-search on the minimum spanning tree from starting node rand_root
        _,parents = sp.sparse.csgraph.breadth_first_order(minT,rand_root, directed =False, return_predecessors =True) 
        
        c[clust_vox] = parents[clust_vox] 
    
    roots = np.where(c==-9999)[0] #we look for -9999 as the root has parents value of -9999, instead of 0.
    c[roots] = roots #unclear what this is doing...replaces the locations of -9999 with the index of where the root is found

    return c



def ComputeCachedLikelihoodTerms(kappa, nu, sigsq):
    cached = [0,kappa, nu, sigsq, nu * sigsq, -gammaln(nu/2) + (1/2)*log(kappa) + (nu/2)*log(nu*sigsq)]
    return cached #of type list


def FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym):  #also implemented in ddCRP
    if sym:
        stats = np.zeros([len(parcels)*(len(parcels)+1)/2,3])
        j = 0 #same as j =1 ?...
        for c1 in range(len(parcels)):
            for c2 in range(c1,len(parcels)):
                #unclear what datatype is D...parcels is a list of list type
                samples = D[parcels[c1],parcels[c2]]
                if c1==c2:
                    # following same as samples = samples(logical(triu(ones(size(samples)),1))); multiplying samples with np.triu(np.ones([samples.shape, samples.shape]),1) essentially zeros out every element not within the triu matrix
                    samples = samples*np.triu(np.ones(samples.shape),1) 
                    if not samples.tolist(): #isempty(samples)
                        continue
                    else:
                        samples = samples.flatten() #samples = samples(:)
                
                stats[j,0] = len(samples)
                stats[j,1] = np.sum(samples)/stats[j,0] #sum(samples)/stats(j,1)
                stats[j,2] = np.sum((samples-stats[j,1])**2) #sum((samples-stats(j,2)).^2)
                j=j+1
        c_sum= sum([1 for i in range(len(c)) if i==c[i]]) #should be same as np.sum(c == 1:len(c)), ie compares 0,1,2,3,4...with vector c and gives 1 if elements are same, then takes sum
        logp = mt.log(alpha) * c_sum + LogLikelihood(stats, hyp)
    else:
        stats = np.zeros([len(parcels)*len(parcels),3])
        j = 0 #same as j = 1?...
        for c1 in range(len(parcels)):
            for c2 in range(len(parcels)):
                samples = D[parcels[c1],parcels[c2]] #unclear of D datatype
                if c1 = c2:                    
                    off_diags = np.ones(samples.shape) #same as true(size(samples))
                    off_diags[0::(samples.shape[0]+1)] = 0 #same as off_diags(1:(size(samples,1)+1):end) = false
                    samples = samples[off_diags]
                    if not samples.tolist():
                        continue
                else:
                    samples = samples.flatten()
                    
                stats[j,0] = length(samples)
                stats[j,1] = np.sum(samples)/stats[j,0] #sum(samples)/stats(j,1)
                stats[j,2] = np.sum((samples-stats[j,1])**2) #sum((samples-stats(j,2)).^2)
                j = j+1
        
        c_sum= sum([1 for i in range(len(c)) if i==c[i]]) #should be same as np.sum(c == 1:len(c)), ie compares 0,1,2,3,4...with vector c and gives 1 if elements are same, then takes sum
        logp = mt.log(alpha) * c_sum + LogLikelihood(stats, hyp) #same as log(alpha) * sum(c' == 1:length(c)) + LogLikelihood(stats, hyp)

    return curr_lp