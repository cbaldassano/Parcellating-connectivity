import numpy as np
import scipy as sp

def LocalSimilarity(D, adj_list, n_clust, vox_to_clust = np.arange(1, D.shape[0]+1)):
       
    N = D.shape[0]
    infn = float("inf")
    W = infn*np.ones((D.shape,D.shape)) #same as inf(size(D))
    
    if CheckSymApprox(D):
        for i in range(D.shape[0]): 
            for j in adj_list[i]:
                #the following lines create same as (1:N ~= i) & (1:N ~= j))
                ij_vector = [1]*N   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
                ij_vector[i]=0
                ij_vector[j]=0  

                jk=(np.arange(N)*ij_vector).tolist()  #this does something like [1,2,3,4,5] * [1,0,1,0,1], and then makes list (list because list is iterable object unlike a numpy array)
                select_vector=[x for x in jk if x!=0]  #this list comprehension only takes the values that are not equal to zero, ie [1,0,3,0,5] becomes [1,3,5]             
                W[i,j] = np.linalg.norm(D[i,select_vector]-D[j,select_vector],2) #2-norm
        for i in range(D.shape[0]): 
            for j in adj_list[i]:
                #the following lines create same as (1:N ~= i) & (1:N ~= j))
                ij_vector = [1]*N   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
                ij_vector[i]=0
                ij_vector[j]=0  

                jk=(np.arange(N)*ij_vector).tolist()  #this does something like [1,2,3,4,5] * [1,0,1,0,1], and then makes list (list because list is iterable object unlike a numpy array)
                select_vector=[x for x in jk if x!=0]  #this list comprehension only takes the values that are not equal to zero, ie [1,0,3,0,5] becomes [1,3,5]
                first_mat=np.concatenate(D[i,select_vector],D[select_vector,i].transpose(),axis=1) #[D(i, (1:N ~= i) & (1:N ~= j)) D((1:N ~= i) & (1:N ~= j), i)']
                second_mat=np.concatenate(D[j,select_vector],D[select_vector,j].transpose(),axis=1)                
                W[i,j] = np.linalg.norm(first_mat-second_mat,2) #2-norm
                
        z = ThresholdSimilarity(W(vox_to_clust,vox_to_clust), n_clust);
    
    return z

def ThresholdSimilarity(W, n_clust):
    thresh = np.mean(W[np.isfinite(W)])
    thresh_step = 0.1
    last_increase = 0
    last_decrease = 0
    
    while(True):
        #following lines is same as G = sparse(W <= thresh) 
        new_W = np.asarray(W)
        low_valued_indices = new_W <= thresh  # Where values are low
        new_W[low_valued_indices] = 0  # All low values set to 0
        G=sp.sparse.csc_matrix(new_W)
        z = graphconncomp(G) 
        curr_n_clust = len(np.unique(z))
        
        if curr_n_clust == n_clust:
            break
        elif curr_n_clust < n_clust:
            if last_increase:
                thresh_step = thresh_step**2
            thresh = thresh * (1-thresh_step)
            last_decrease = 1
            last_increase = 0
        else:
            if last_decrease:
                thresh_step = thresh_step**2
            thresh = thresh * (1+thresh_step)
            last_decrease = 0
            last_increase = 1
            
    return z
    
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
    
def graphconncomp(G):
    _, z = sp.sparse.csgraph.connected_components(G,directed=False,connection='weak',return_labels=True)
    return z