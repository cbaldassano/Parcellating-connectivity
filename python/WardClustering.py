import numpy as np
import scipy as sp

def WardClustering(D, adj_list, n_clust, vox_to_clust):
    if nargin == 4: #no nargin in python
        size = D.shape[0] #this is equivalent to size(D,1)
        vector = np.arange(size) #this creates a vector that is equivalent to 1:size(D,1)
        
        #for next line, what data type is vox_to_clust...assuming it is an array. Converts vector and vox_to_list to list first to do comparision            
        
        unclust=[1*(not(vector.tolist() == vox_to_clust.tolist()) for x in range(len(A))]  #this creates list comprehension if equal, then it is 0, else 1
        unclust = np.array(unclust)    #convert unclust to type array            
        adj_list = []
        n_clust = n_clust + np.sum(unclust)
        
    Z = LinkageConstrained(D, adj_list)
    z= sp.cluster.hierarchy.fcluster(Z, t=n_clust, criterion = 'maxclust')  #should be same as z = cluster(Z, 'maxclust', n_clust)    
    
    return(z,Z)
