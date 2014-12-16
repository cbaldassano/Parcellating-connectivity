import numpy as np
import random

def GenerateRandomClustering(adj_list, K, z_init):
    N = len(adj_list)
    if nargin < 3: #no nargin in python
         merge_vox = np.arange(1,N,1) #vector 1 through N
         K_curr = N
         z = np.arange(1,N,1) #vector 1 through N
    else:
         z=z_init
         merge_vox = np.where(z_init==0)   #equivalent to find in matlab
         K_curr = len(merge_vox);
         z[merge_vox] = np.arange((z_init.max(0)+1),(z_init.max(0)+K_curr),1)
         
    while(K_curr > K):
         clustLabels = np.unique(z[merge_vox])
         n=random.randint(1, len(clustLabels))  #random number over len of clustlabels
         zToMerge = clustLabels(np.array([[n]]))  #why creating a 1x1 random matrix?  over len of clustlabels
         #not sure what this does ??mergeChoices = unique(z([adj_list{z==zToMerge}]))
         mergeChoices= np.array(list(set(mergeChoices) - set(zToMerge)))
         if mergeChoices == True:
             n=random.randint(1, len(mergeChoices))  #random number over len of clustlabels
             z(z==zToMerge) = mergeChoices(np.array([[n]]))   #Unclear what z(z==zToMerge) represents
             K_curr = K_curr - 1
             
    zRelabel = np.zeros((3,4))
    relabelInd = 1
    
    for v in range(np.unique(z)):
        zRelabel(z==v) = relabelInd
        relabelInd = relabelInd+1
        
    z = zRelabel
    if nargout == 2:   #not clear here of how to determine number of output
        c = ClusterSpanningTrees(z, adj_list)
    return(z,c)