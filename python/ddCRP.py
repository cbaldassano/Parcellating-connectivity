import numpy as np
from scipy import special
from scipy import sparse
import random as rd
import math as mt
import time
import StatsUtil


from pylab import *
def PlotC(c):
    for i in range(len(c)):
        x1 = mt.floor(i/18)
        y1 = i-x1*18
        x2 = mt.floor(c[i]/18)
        y2 = c[i]-x2*18

        plot([x1, x2], [y1, y2], color='k', linestyle='-', linewidth=2)
def PlotAdj(adj_list):
    for i in range(len(adj_list)):
        x1 = mt.floor(i/18)
        y1 = i-x1*18
        for j in adj_list[i]:
            x2 = mt.floor(j/18)
            y2 = j-x2*18
            plot([x1, x2], [y1, y2], color='k', linestyle='-', linewidth=2)


def ddCRP(D, adj_list, init_c, gt_z, num_passes, alpha, kappa, nu, sigsq, stats_interval, verbose):
    map_z =  np.zeros(np.shape(D)[0])
    stats = {'times':[], 'lp':[], 'NMI':[], 'K':[], 'z':[], 'c':[]}
    
    hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq)
    nvox=len(adj_list)
    
    # Generate random initialization if not specified
    if init_c.size==0:
        c = np.zeros(nvox)
        for i in range(nvox):
            neighbors = np.concatenate((adj_list[i], i),axis=1) 
            c[i] = neighbors[rd.randint(1,len(neighbors))]
    else:
        c = init_c
            
    G = sparse.coo_matrix((np.ones(nvox),(np.arange(nvox),c)), shape=(nvox,nvox))
    K, z, parcels = ConnectedComp(G)
    
    sym = StatsUtil.CheckSymApprox(D)
    curr_lp = FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym)
    
    max_lp = -float('inf')
    steps = 0
    t0 = time.clock()
    
    for curr_pass in range(num_passes):
        order = np.random.permutation(nvox)
        
        for i in order:
            if curr_lp > max_lp:
                max_lp = curr_lp
                map_z = z
            
            if steps % stats_interval == 0:
                stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, verbose);
        
            # Compute change in log-prob when removing the edge c_i
            CooModifyRow(G, i, -1)
            if c[i] == i:
                rem_delta_lp, z_rem, parcels_rem = -mt.log(alpha), z, parcels
            else:
                K_rem, z_rem, parcels_rem = ConnectedComp(G)
                
                if K_rem != K:
                    rem_delta_lp = -LikelihoodDiff(D, parcels_rem, z_rem[i], z_rem[c[i]], hyp, sym)
                else:
                    rem_delta_lp = 0

            # Compute change in log-prob for each possible edge c_i
            adj_list_i = adj_list[i]
            lp = np.zeros((len(adj_list_i)+1))
            lp[len(adj_list_i)] = mt.log(alpha)
            cached_merge = -1*np.ones(len(adj_list_i), dtype=int32)
            for n_ind in range(len(adj_list_i)):
                n = adj_list_i[n_ind]
                if z_rem[n] == z_rem[c[i]]:
                    # Just undoing edge removal
                    lp[n_ind] = -rem_delta_lp - (c[i] == i)*mt.log(alpha)
                elif z_rem[n] != z_rem[i]:
                    # Proposing novel merge
                    # First check cache to see if this is already computed
                    prev_lp = np.flatnonzero(cached_merge == z_rem[n])
                    if prev_lp.size > 0:
                        lp[n_ind] = lp[prev_lp[0]]
                    else:
                        lp[n_ind] = LikelihoodDiff(D, parcels_rem, z_rem[i], z_rem[n], hyp, sym)
                        cached_merge[n_ind] = z_rem[n]
            
            # Pick new edge proportional to probability
            new_neighbor = ChooseFromLP(lp)
            if new_neighbor < len(adj_list_i):
                c[i] = adj_list_i[new_neighbor]
            else:
                c[i] = i
                
            curr_lp = curr_lp + rem_delta_lp + lp[new_neighbor]
            CooModifyRow(G, i, c[i])
            K, z, parcels = ConnectedComp(G)
            steps = steps + 1

            
    stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, verbose)
    return (map_z, stats)


def ConnectedComp(G):
    # Compute connected components (number and component labels)
    K, z = sparse.csgraph.connected_components(G,directed=False,connection='weak',return_labels=True)

    sorted_i = np.argsort(z)
    sorted_z = np.sort(z)
    parcels = np.split(sorted_i,np.flatnonzero(np.diff(sorted_z))+1)
        
    return(K, z, parcels)
    

def LikelihoodDiff(D, parcels_split, split_i1, split_i2, hyp, sym):
    
    # Compute log-likelihood for split parcels
    K = len(parcels_split)
    s=np.zeros((K,K,3))
    for split_ind in [split_i1, split_i2]:

        # Get all connections to split_ind parcel
        for i in range(K):
            samples = D[np.ix_(parcels_split[i],parcels_split[split_ind])]
            if i == split_ind:
                if sym:
                    samples = samples[np.triu_indices(len(samples),1)]
                else:
                    samples = [];
            else: 
                samples = samples.ravel()
            s[i,split_ind,:] = SufficientStats(samples)  
        
        # If not symmetric, get all connections from split_ind parcel
        if not sym:
            for i in range(K):
                samples = D[np.ix_(parcels_split[split_ind], parcels_split[i])]
                if i == split_ind:
                    off_diags = np.logical_not(np.eye(np.shape(samples)[0],dtype='bool'))
                    samples = samples[off_diags]
                else:
                    samples = samples.ravel()
                s[split_ind,i,:] = SufficientStats(samples)
    
    # Compute log-likelihood of all sufficient stats for split parcels
    if sym:
        suffstats1 = np.reshape(s[:,split_i1,:],(-1,3))
        suffstats2 = np.delete(np.reshape(s[:,split_i2,:],(-1,3)),split_i1,0)
        split_ll = LogLikelihood(np.concatenate((suffstats1,suffstats2)), hyp)
    else:
        suffstats_to1 = np.reshape(s[:,split_i1,:],(-1,3))
        suffstats_to2 = np.reshape(s[:,split_i2,:],(-1,3))
        suffstats_from1 = np.delete(np.reshape(s[split_i1,:,:],(-1,3)),[split_i1, split_i2],0)
        suffstats_from2 = np.delete(np.reshape(s[split_i2,:,:],(-1,3)),[split_i1, split_i2],0)
        split_ll = LogLikelihood(np.concatenate((suffstats_to1, suffstats_to2, suffstats_from1, suffstats_from2)), hyp)    
    

    # Compute log-likelihood for all merged parcels
    m = np.zeros((2, K, 3))
    # Compute merges for all off-diagonal parcels
    for dir in range(2):
        if dir ==1:
            if sym:
                break
            else:
                s = np.transpose(s, (1, 0, 2))

        for i in range(K):
            if i != split_i1 and i != split_i2:
                s_m = np.reshape(s[i, (split_i1, split_i2), :], (2,3))
                m[dir,i,:] = MergeSuffStats(s_m)
    if sym:
        # Compute central diagonal merge
        m_central = MergeSuffStats(np.concatenate((
            MergeSuffStats(np.reshape(s[split_i1, (split_i1, split_i2), :], (2,3))),
            np.reshape(s[split_i2, split_i2, :], (1,3)))))

        # Compute log-likelihood of all sufficient stats for merged parcels
        merge_ll = LogLikelihood(np.concatenate((
            np.reshape(m[0,:,:], (-1,3)),
            m_central)), hyp)
    else:
        # Compute central diagonal merge
        m_central = MergeSuffStats(np.concatenate((
            MergeSuffStats(np.reshape(s[split_i1, (split_i1, split_i2), :], (2,3))).reshape((1,3)),
            MergeSuffStats(np.reshape(s[split_i2, (split_i1, split_i2), :], (2,3))).reshape((1,3)))))

        # Compute log-likelihood of all sufficient stats for merged parcels
        merge_ll = LogLikelihood(np.concatenate((
            np.reshape(m[0,:,:], (-1,3)),
            np.reshape(m[1,:,:], (-1,3)),
            m_central.reshape((1,3)))), hyp)
    
    ld = merge_ll - split_ll                    
    return ld

# Compute count, mean, and num of squared deviations
def SufficientStats(samples):
    suffstats = np.zeros(3)
    
    if len(samples)==0: 
        return suffstats
    
    suffstats[0] = len(samples)
    suffstats[1] = np.sum(samples)/suffstats[0]
    suffstats[2] = np.sum((samples-suffstats[1])**2)
    return suffstats

# Compute sufficient statistics for merging two sets of suff stats
def MergeSuffStats(s_m):
    m = np.zeros(3)
    m[0] = s_m[0,0] + s_m[1,0]
    m[1] = (s_m[0,0]*s_m[0,1] + s_m[1,0]*s_m[1,1])/m[0]
    m[2] = s_m[0,2] + s_m[1,2] + (s_m[0,0]*s_m[1,0])/m[0] * (s_m[0,1] - s_m[1,1])**2
    return m
    
# Update diagnostic stats
def UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, verbose):    
    stats['lp'].append(curr_lp)
    stats['K'].append(K)
    stats['z'].append(z)
    stats['c'].append(c)
    curr_time = time.clock() - t0
    stats['times'].append(curr_time)
    if verbose:
        print('Step: ' + str(steps) + ' Time: ' + str(curr_time) + ' LP: ' + str(curr_lp) + ' K: ' + str(K))

    if gt_z.size > 0:
        stats['NMI'].append(StatsUtil.NMI(gt_z, map_z))

    return stats

# Precompute and package hyperparameter expressions into hyp
def ComputeCachedLikelihoodTerms(kappa, nu, sigsq):
    cached = [0,kappa, nu, sigsq, nu * sigsq, -mt.lgamma(nu/2) + (1/2)*mt.log(kappa) + (nu/2)*mt.log(nu*sigsq)]
    return cached


# Compute full probability of a given parcellation
# Slow, should only be used during initialization (updated incrementally during inference)
def FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym):
    self_loops = sum([1 for i in range(len(c)) if i==c[i]])

    if sym:
        stats = np.zeros([len(parcels)*(len(parcels)+1)/2,3])
        j = 0
        for c1 in range(len(parcels)):
            for c2 in range(c1,len(parcels)):
                samples = D[np.ix_(parcels[c1],parcels[c2])]
                if c1==c2:
                    samples = samples[np.triu_indices(len(samples),1)]
                    if len(samples)==0:
                        continue
                else:
                    samples = samples.ravel()
                
                stats[j,0] = len(samples)
                stats[j,1] = np.sum(samples)/stats[j,0]
                stats[j,2] = np.sum((samples-stats[j,1])**2)
                j=j+1
        
        logp = mt.log(alpha) * self_loops + LogLikelihood(stats, hyp)
    else:
        stats = np.zeros([len(parcels)*len(parcels),3])
        j = 0
        for c1 in range(len(parcels)):
            for c2 in range(len(parcels)):
                samples = D[np.ix_(parcels[c1],parcels[c2])]
                if c1 == c2:
                    off_diags = np.logical_not(np.eye(np.shape(samples)[0],dtype='bool'))
                    samples = samples[off_diags]
                    if len(samples)==0:
                        continue
                else:
                    samples = samples.ravel()
                    
                stats[j,0] = len(samples)
                stats[j,1] = np.sum(samples)/stats[j,0]
                stats[j,2] = np.sum((samples-stats[j,1])**2)
                j = j+1
        
        logp = mt.log(alpha) * self_loops + LogLikelihood(stats, hyp)

    return logp

# Compute log-likelihood of model given suff stats
def LogLikelihood(stats, hyp):
    stats = stats[stats[:,0]>1,:]
    if stats.size > 0:
        # stats = [N | mu | sumsq]
        # hyp = [mu0 kappa0 nu0 sigsq0 nu0*sigsq0 const_logp_terms] as type list
        kappa = hyp[1] + stats[:,0]
        nu = hyp[2] + stats[:,0]
        # Assume mu0=0 and kappa0 << n
        nu_sigsq = hyp[4] + stats[:,2] + hyp[1] * (stats[:,1])**2
        
        #logp = sum(hyp(6) + gammaln(nu/2)- 0.5*log(kappa) - (nu/2).*log(nu_sigsq)- (stats(:,1)/2)*log(pi))    
        logp = np.sum(hyp[5] + special.gammaln(nu/2)- 0.5*np.log(kappa) - \
            (nu/2)*np.log(nu_sigsq)- (stats[:,0]/2)*np.log(mt.pi))
  
        return logp
    else:
        return 0

# Randomly choose index based on unnormalized log probabilities
def ChooseFromLP(lp):
    import random
    max_lp = lp.max()
    normLogp = lp - (max_lp + np.log(np.sum(np.exp(lp-max_lp))))
    p = np.exp(normLogp)
    p[np.isfinite(p)==False]=0
    cumP = np.cumsum(p)
    i = np.where(cumP>random.random())[0][0]
    return i

# Insert value in coo_matrix at row,col, or remove row if col==-1
def CooModifyRow(G,row,col):
    if col == -1:
        G.row = np.delete(G.row,row)
        G.col = np.delete(G.col,row)
        G.data = np.delete(G.data,row)
    else:
        G.row = np.insert(G.row, row, row)
        G.col = np.insert(G.col, row, col)
        G.data = np.insert(G.data, row, 1)
    return G