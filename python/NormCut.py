import numpy as np

def NormCut(D, adj_list, n_clust, vox_to_clust = np.arange(D.shape[0])):
  
    # not implemented : addpath('normcut')
    N = D.shape[0]
    
    W = np.zeros((D.shape, D.shape))
    eps = 0.0001
    if CheckSymApprox(D):
        for i in range(N):
            for j in adj_list[i]:
                #the following lines create same as (1:N ~= i) & (1:N ~= j))
                ij_vector = [1]*N   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
                ij_vector[i]=0
                ij_vector[j]=0  

                jk=(np.arange(N)*ij_vector).tolist()  #this does something like [1,2,3,4,5] * [1,0,1,0,1], and then makes list (list because list is iterable object unlike a numpy array)
                select_vector=[x for x in jk if x!=0]  #this list comprehension only takes the values that are not equal to zero, ie [1,0,3,0,5] becomes [1,3,5]             
                W[i,j] = 1/((np.linalg.norm(D[i,select_vector]-D[j,select_vector],2))+eps)
                
        for i in range(N):
            for j in adj_list[i]:
                #the following lines create same as (1:N ~= i) & (1:N ~= j))
                ij_vector = [1]*N   
                ij_vector[i]=0
                ij_vector[j]=0  

                jk=(np.arange(N)*ij_vector).tolist() 
                select_vector=[x for x in jk if x!=0]  
                first_mat=np.concatenate(D[i,select_vector],D[select_vector,i].transpose(),axis=1) #[D(i, (1:N ~= i) & (1:N ~= j)) D((1:N ~= i) & (1:N ~= j), i)']
                second_mat=np.concatenate(D[j,select_vector],D[select_vector,j].transpose(),axis=1)                
                W[i,j] = 1/((np.linalg.norm(first_mat-second_mat,2))+eps)
                
    nc_discrete = ncutW(W[vox_to_clust,vox_to_clust], n_clust) #ncutW will have to be implemented.  Looks as large library written by some other group.  Maybe shortcut?
    z = nc_discrete * np.arange(n_clust) #is n_clust a scalar? transpose not taken as python defaults as column vectors
    return z
                