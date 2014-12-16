import numpy as np
import scipy as sp
import collections
import random as rd
import math as mt

def ddCRP(D, adj_list, init_c, gt_z, num_passes, alpha, kappa, nu, sigsq, stats_interval, verbose):
    map_z =  np.zeros(np.shape(D)[0])
    StatStruct = collections.namedtuple('Stats',['times','lp','NMI','K','z','c'])
    stats = StatStruct
    
    
    hyp = ComputeCachedLikelihoodTerms(kappa, nu, sigsq) #hyp as type list
    pid = rd.randint(1,10000)
    nvox=len(adj_list)
    
    #what data type is const_c, init_c?
    if not const_c: const_c = np.zeros(nvox)
    if not init_c: init_c = np.zeros(nvox)
    if burn_in_passes: pair_prob = np.zeros(nvox*(nvox-1)/2):
    
    c = const_c
    
    for i in np.where(const_c==0)[0]: #where returns tuple
        if init_c[i]==0:
            #is neighbors an array or other data structure?
            neighbors = np.concatenate((adj_list[i], i),axis=1) 
            c[i] = neighbors[rd.randint(1,len(neighbors))]
        else:
            c[i] = init_c[i]
            
    G = sp.sparse.csc_matrix((np.ones(nvox),(np.arange(nvox),c)), shape=(nvox,nvox)) #same as G = sparse(1:nvox,c,1,nvox,nvox). Must use np.ones(nvox) and not 1
    K, z, parcels = ConnectedComp(G)
    
    sym = CheckSymApprox(D)
    curr_lp = FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym)
    
    max_lp = -float('inf')
    #tic toc not implemented yet 
    steps = 0
    
    for pass1 in range(num_passes):  #"pass" changed to "pass1" as "pass" is reserved word in python
        nonconst_vox = np.where(const_c==0)[0]
        order = nonconst_vox[np.random.permutation(len(nonconst_vox))] 
        
        #what datatype is order?  
        for i in order:
            if curr_lp > max_lp:
                max_lp = curr_lp
                map_z = z
            
            if burn_in_passes and pass1 > burn_in_passes: # if (~isempty(burn_in_passes) && pass > burn_in_passes)
                pair_prob = pair_prob + (1 - sp.spatial.distance.pdist(z, 'hamming')) 
            if mod[steps,stats_interval] == 0:
                stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, pid, verbose);
        
            if c[i] == i:
                rem_delta_lp = -mt.log(alpha)
                z_rem = z 
                parcels_rem = parcels
            else:
                G[i,c[i]] = 0
                K_rem, z_rem, parcels_rem = ConnectedComp(G)
                
                if K_rem != K:
                    rem_delta_lp = -LikelihoodDiff(D, parcels_rem, z_rem[i], z_rem[c[i]], hyp, sym)
                else:
                    rem_delta_lp = 0
            
            adj_list_i = adj_list[i]
            lp = np.zeros((len(adj_list_i)+1))
            lp[len(adj_list_i)] = mt.log(alpha) #end should be equal to len(adj_list_i)
            
            for n_ind in range(len(adj_list_i)):
                n = adj_list_i[n_ind]
                if z_rem[n] == z_rem[c[i]]: #Clustered with old neighbor
                    lp[n_ind] = -rem_delta_lp
                elif z_rem[n] != z_rem[i]:  #Not already clustered with n
                    lp[n_ind] = LikelihoodDiff(D, parcels_rem, z_rem[i], z_rem[n], hyp, sym)
            
            new_neighbor = ChooseFromLP(lp)
            if new_neighbor <= len(adj_list_i):
                c[i] = adj_list_i[new_neighbor]
            else:
                c[i] = i
                
            curr_lp = curr_lp + rem_delta_lp + lp[new_neighbor]
            G[i,c[i]] = 1
            K, z, parcels = ConnectedComp(G)
            steps = steps + 1
            
    stats = UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, pid, verbose)    
    
    if burn_in_passes: # ie ~isempty(burn_in_passes):
        pair_prob = pair_prob / (np.sum((const_c==0)*(num_passes-burn_in_passes)))
    else:
        pair_prob = np.array([])
        
    return (map_z, stats, pair_prob)

def ConnectedComp(G): 
    K, z = sp.sparse.csgraph.connected_components(G,directed=False,connection='weak',return_labels=True)

    sorted_z = np.sort(z, axis=0)
    sorted_i = np.argsort(z, axis=0) #the indices when sorted
         
    #following lines creates [0 sorted_z (K+1)]
    diff_input = np.append(sorted_z,K+1)
    diff_input = np.insert(diff_input,0,0)
    
    find_diff = np.nonzero(np.diff(diff_input)) #find(diff([0 sorted_z (K+1)]))
    diff1 = np.diff(find_diff) #diff(find(diff([0 sorted_z (K+1)])))
    
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
        
    return(K, z, parcels)

def CheckSymApprox(D):
    #equivalent to sym_sub = [randi(size(D,1), 1000,1) randi(size(D,1), 1000,1)
    #sym_sub is a 1000 x 2 random integer matrix
    sym_sub = np.random.randint(D.shape[0], size=(1000,2)) 
    
    #following should be same as sub2ind(size(D), sym_sub(:,1), sym_sub(:,2)) and sub2ind(size(D), sym_sub(:,2), sym_sub(:,1))
    a = np.ravel_multi_index((sym_sub[:,0],sym_sub[:,1]), dims=np.shape(D), order='F') #'F' represents Fortran order indexing which matlab is based upon.
    b = np.ravel_multi_index((sym_sub[:,1],sym_sub[:,0]), dims=np.shape(D), order='F')
    
    #following is sym = all(D(...)==D(...)
    sym = np.all(D.ravel()[a] == D.ravel()[b])
    
    return sym
    
def LikelihoodDiff(D, parcels_split, split_i1, split_i2, hyp, sym):
    K = len(parcels_split)
    s=np.zeros((3,K,K)) #different than matlab which is s = zeros(K, K, 3)
    #ASSUMPTION: split_i1 and split_i2 are vectors
    split_ind1  = np.concatenate((split_i1, split_i2), axis=1)    
    split_ind_list = split_ind1.tolist() #change split_ind from vector to list so that it is iterable in for loop    
    
    for split_ind in split_ind_list: #split_ind is element of list split_ind_list
        for i in range(K):
            samples = D[parcels_split[i], parcels_split[split_ind]]
            if sym and i == split_ind:
                b=np.triu(np.ones(len(samples),len(samples)))>0 #creates matrix of type boolean. b is matrix of type bool
                samples = samples(b) #samples is now a vector of all nonzero values
            else: 
                samples = samples.flatten() #this rolls out the matrix into 1-dim vector. same as  samples = samples(:)
            s[:,i,split_ind] = SufficientStats(samples)  
            
        if not sym:
            for i in range(K):
                samples = D[parcels_split[split_ind], parcels_split[i]]
                if i == split_ind:
                    off_diags = np.ones(samples.shape) #in python, shape gives same information as size. Matrix of ones is same as matlab's True matrix
                    off_diags[0::(samples.shape[0]+1)] = 0 #same as off_diags(1:(size(samples,1)+1):end) = false
                    samples = samples[off_diags]
                else:
                    samples = samples.flatten()
                s[:,split_ind,i] = SufficientStats(samples)
                
    if sym:
        #Probably wrong...Unclear what is meant by K ~= split_i1, split_i2
        #[] from matlab represents size of the 2-Dim matrix of 3-Dim s matrix.   
        s_dim = (s.size)/3    #[] should equal s_dim??...
        split_ll = LogLikelihood(np.concatenate([np.reshape(s([::,split_i1],s_dim,3), np.reshape(s[0::K != split_i1 and K !=split_i2],s_dim,3)],axis=0), hyp)
    else:
        #following lines breaks down split_ll = LogLikelihood([reshape(s(:, split_i1,:),[],3); reshape(s(:, split_i2,:),[],3); reshape(s(split_i1, (1:K ~= split_i1) & (1:K ~= split_i2),:),[],3); reshape(s(split_i2, (1:K ~= split_i1) & (1:K ~= split_i2),:),[],3)], hyp)
        first_concat = np.concatenate([np.reshape(s[:, split_i1,:],s_dim,3), np.reshape(s[:, split_i2,:],s_dim,3)], axis=0)
        
        #the following lines create same as (1:K ~= split_i1) & (1:K ~= split_i2))
        ij_vector = [1]*K   #created list as opposed to using ones(K) vector because multiplication with np.arange gives a dimension for list that can't be iterated
        ij_vector[split_i1]=0
        ij_vector[split_i2]=0    
            
        second_concat = np.concatenate([np.reshape(s[:,split_i1, ij_vector],s_dim,3), np.reshape(s[:,split_i2, ij_vector],s_dim,3)], axis=0)
        final_concat = np.concatenate([first_concat, second_concat], axis = 0) #final concatenation
        
        split_ll  = LogLikelihood(final_concat,hyp)
    
    m = np.zeros(3, 2, K) #same as zeros(2, K, 3)
    
    for dir in range(2)
        if dir ==1: #same as dir == 2 in matlab
            if sym:
                break
            else:
                np.transpose(2,(0,2,1)) #Same as s = permute(s, [2 1 3]).  Assuming you are just taking transpose in 2 dimensions
        for i in range(K):
            if i != split_i1 and i != split_i2:
                split_list = np.concatenate([split_i1, split_i2], axis = 1)
                split_list = split_list.tolist()
                s_m = np.reshape(s[:, i, split_list],2,3) #split_list is a list that selects the columns of [split_i1 split_i2]
                m[:,dir,i] = MergeSuffStats(s_m)
    if sym:
        #following same as m_central = MergeSuffStats(...
        #..[MergeSuffStats(reshape(s(split_i1, [split_i1 split_i2], :),2,3)); 
        #..reshape(s(split_i2, split_i2, :),1,3)])
        split_list = np.concatenate([split_i1, split_i2], axis = 1)
        split_list = split_list.tolist() #same as [split_i1 split_i2], except list
        a = MergeSuffStats(np.reshape(s[:, split_i1.tolist(), split_list],2,3)) #same as MergeSuffStats(reshape(s(split_i1, [split_i1 split_i2], :),2,3))
        b = np.reshape(s[:,split_i2.tolist(), split_i2.tolist()],1,3)  #same as reshape(s(split_i2, split_i2, :),1,3)
        merge_input = np.concatenate([a,b], axis = 0) #same as [MergeSuffStats(reshape(s(split_i1, [split_i1 split_i2], :),2,3)); reshape(s(split_i2, split_i2, :),1,3)]
        m_central = MergeSuffStats(merge_input)
        
        c =np.reshape(m[:,0,:], (m.size)/3, 3) #Not sure if [] equals (m.size)/3 here
        d = np.concatenate([c,m_central], axis = 0) # d represents the [reshape(m(1, :, :),[],3); m_central]
        merge_ll = LogLikelihood(d, hyp)
    else:
        #following same as  m_central = MergeSuffStats([MergeSuffStats(reshape(s(split_i1, [split_i1 split_i2], :),2,3)); 
        #...MergeSuffStats(reshape(s(split_i2, [split_i1 split_i2], :),2,3))]);
        split_list = np.concatenate([split_i1, split_i2], axis = 1)
        split_list = split_list.tolist() #same as [split_i1 split_i2], except list
        a = np.reshape(s[:,split_i1, split_list],2,3) #same as reshape(s(split_i1, [split_i1 split_i2], :),2,3)
        b= np.reshape(s[:,split_i2, split_list],2,3) #same as reshape(s(split_i2, [split_i1 split_i2], :),2,3))
        merge_input = np.concatenate([MergeSuffStats(a),MergeSuffStats(b)], axis = 0)
        m_central = MergeSuffStats(merge_input)
        
        #following same as merge_ll = LogLikelihood([reshape(m(1, :, :),[],3); reshape(m(2, :, :),[],3)
        a = reshape(m[:, 0, :],(m.size)/3, 3)
        b = reshape(m[:, 1, :],(m.size)/3, 3)
        log_input = np.concatenate([a,b, m_central], axis = 0) #same as [reshape(m(1, :, :),[],3); reshape(m(2, :, :),[],3); m_central]
        merge_ll = LogLikelihood(log_input, hyp)
    
    ld = merge_ll - split_ll                    
    return ld

def SufficientStats(samples):
    suffstats = np.zeros(3)
    #is samples a vector?
    samples_list = samples.tolist()
    
    if not samples_list: 
        return
    
    suffstats[0] = len(samples)
    suffstats[1] = np.sum(samples)/suffstats[0]
    suffstats[2] = np.sum((samples-suffstats[1])**2)
    return suffstats

def MergeSuffStats(s_m):
    m = np.zeros(3)
    m[0] = s_m[0,0] + s_m[1,0] #same as m(1) = s_m(1,1) + s_m(2,1)
    m[1] = (s_m[0,0]*s_m[0,1] + s_m[1,0]*s_m[1,1])/m[0] #same as  m(2) = (s_m(1,1)*s_m(1,2) + s_m(2,1)*s_m(2,2))/m(1)
    m[2] = s_m[0,2] + s_m[1,2] + (s_m[0,0]*s_m[1,0])/m[0] * (s_m[0,1] - s_m[1,1])**2 #same as m(3) = s_m(1,3) + s_m(2,3) + (s_m(1,1)*s_m(2,1))/m(1) * (s_m(1,2) - s_m(2,2))^2    
    return m
    
def UpdateStats(stats, t0, curr_lp, K, z, c, steps, gt_z, map_z, pid, verbose):    
    stats.lp.append(curr_lp) #stats.lp = [stats.lp curr_lp]
    stats.K.append(K) #stats.K = [stats.K K]
    stats.z.append(z) #stats.z = [stats.z; z]
    #elapsed = toc(t0) not implemented
    #stats.times = [stats.times elapsed]  not implemented
    stats.c.append(c)
    #following not implemented
    #if (verbose)
    #    disp(['Step: ' num2str(steps) ...
    #          '  Time: ' num2str(elapsed) ...
    #          '  LP: ' num2str(curr_lp) ...
    #          '  K: ' num2str(K)]);
    #end
    #if (~isempty(gt_z))
    #    stats.NMI = [stats.NMI CalcNMI(gt_z, map_z)];
    #end
    #save(['/data/supervoxel/output/temp/' num2str(pid) '.mat'], ...
    #    'map_z', 'stats');
    return stats
    
def ComputeCachedLikelihoodTerms(kappa, nu, sigsq):
    cached = [0,kappa, nu, sigsq, nu * sigsq, -gammaln(nu/2) + (1/2)*log(kappa) + (nu/2)*log(nu*sigsq)]
    return cached #of type list


def FullProbabilityddCRP(D, c, parcels, alpha, hyp, sym):
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
        c_sum= sum([1 for i in range(len(c)) if i==c[i]]) #should be same as sum(c' == 1:length(c)), ie compares 0,1,2,3,4...with vector c and gives 1 if elements are same, then takes sum
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
    
def LogLikelihood(stats, hyp):
    stats = stats[stats[:,0]>1,:] #selects the nonzero rows of stats; selects those rows of stats for which stats[:,0]>1...the first column elements are greater than 1
    # stats = [N | mu | sumsq]
    # hyp = [mu0 kappa0 nu0 sigsq0 nu0*sigsq0 const_logp_terms] as type list
    kappa = hyp[1] + stats[:,0] #kappa = hyp(2) + stats(:,1)
    nu = hyp[2] + stats[:,0] #hyp(3) + stats(:,1)
    # nu_sigsq = hyp(5) + sumSqX + (n*hyp(2)) / (hyp(2)+n) * (hyp(1) - meanX)^2;
    # Assume mu0=0 and kappa0 << n
    nu_sigsq = hyp[4] + stats[:,2] + hyp[1] * (stats[:,1])**2 #nu_sigsq = hyp(5) + stats(:,3) + hyp(2) * stats(:,2).^2
    
    #logp = sum(hyp(6) + gammaln(nu/2)- 0.5*log(kappa) - (nu/2).*log(nu_sigsq)- (stats(:,1)/2)*log(pi))    
    logp = np.sum(hyp[5] + mt.lgamma(nu/2)- 0.5*mt.log(kappa) - (nu/2)*mt.log(nu_sigsq)- (stats[:,0]/2)*mt.log(mt.pi))
  
    return logp

def ChooseFromLP(lp):
    import random
    from numpy import inf
    max_lp = lp.max(0)
    f=np.vectorize(lambda x: mt.exp(x)) #f is a function of exp(x), except you can overload x with vector, ie elementwise operation on vector
    normLogp = lp - (max_lp + mt.log(np.sum(f(lp-max_lp)))) #python you can't pass a vector argument to a operation like exponentiation
    p = f(normLogp) 
    p[np.isfinite(p)==False]=0 #equivalent to p(~isfinite(p)) = 0
    cumP = np.cumsum(p)
    i = np.where(cumP>random.random)[0][0] #equivalent to i = find(cumP>rand,1), first instance    
    return i
  

if __name__ == "__main__":
    res = ddCRP(np.zeros((10,10)),0,0,0,0,0,0,0,0,0,0)
    print(res[0])
