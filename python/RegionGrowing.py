import numpy as np

def RegionGrowing(subject, experiment, n_clust):
    #will have to implement loading data    
    #stats = struct('NMI',[], 'conn_diff',[]);
    #loaded = load(['../data/' subject '/' experiment '.mat']);
    #D = loaded.D;
    #adj_list = loaded.adj_list;

    N=D.shape[0]
    
    if experiment == "PPA":
        #labels = loaded.labels;
        #bold = loaded.bold;
    
    if subject == "synth":
        #gt_z = loaded.z;
    else:
        gt_z = np.array([]) #datatype? 
    
    #clear loaded
    
    instability = np.zeros(N)
    
    for i in range(N):
        for j in adj_list[i]:
            #the following lines create same as (1:N ~= i) & (1:N ~= j))
            ij_vector = [1]*N   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
            ij_vector[i]=0
            ij_vector[j]=0  

            jk=(np.arange(N)*ij_vector).tolist()  #this does something like [1,2,3,4,5] * [1,0,1,0,1], and then makes list (list because list is iterable object unlike a numpy array)
            select_vector=[x for x in jk if x!=0]  #this list comprehension only takes the values that are not equal to zero, ie [1,0,3,0,5] becomes [1,3,5]
            instability[i] = instability[i] +np.linalg.norm(D[i,select_vector]-D[j,select_vector],2)
        instability[i] = instability[i] / len(adj_list[i])
        
    for i in range(N):
        instability[i] = np.mean([instability[i], np.mean(instability[adj_list[i]])])
    
    seeds = [] #datatype?
    parcels = [] #datatype? list datatype okay?
    parcel_feat = np.zeros(N) #"zeros(0,N)" creates a 0x10 matrix...is this correct?
    parcel_nbr = [] #datatype? list datatype okay?
    parcel_nbr_diff = []
    added_to_parcel = np.zeros(N) 
    
    for i in range(N):
        if np.all(instability(i) <= instability(adj_list[i])):
            seeds.append(i)
            parcels.append(i)
            added_to_parcel[i] = 1
            instability[i] = -1
            parcel_feat= np.column_stack((parcel_feat,D[i,:])).transpose() #for whatever reason concatenate won't work when array is (n,~) in dimension.
            parcel_nbr = parcel_nbr.append(adj_list[i])
            nbr_diff = np.zeros(len(adj_list[i]))
            
            for j in range(len(adj_list[i])):
                #the following lines create same as (1:N ~= i) & (1:N ~= adj_list{i}(j)))
                ij_vector = [1]*N   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
                ij_vector[i]=0
                ij_vector[adj_list[i][j]]=0  

                jk=(np.arange(N)*ij_vector).tolist()  #this does something like [1,2,3,4,5] * [1,0,1,0,1], and then makes list (list because list is iterable object unlike a numpy array)
                select_vector=[x for x in jk if x!=0]  #this list comprehension only takes the values that are not equal to zero, ie [1,0,3,0,5] becomes [1,3,5]             
                nbr_diff[j] = np.linalg.norm(D[i,select_vector]-D[adj_list[i][j],select_vector],2)
            parcel_nbr_diff = parcel_nbr_diff.append(nbr_diff)
    
    while np.sum(added_to_parcel) < N:
        max_diff = np.amax(np.asarray(parcel_nbr_diff))
        min_diff = np.amin(np.asarray(parcel_nbr_diff))
        thresh = min_diff + 0.1*(max_diff - min_diff)
        
        add_mask = [[[0]*shape(parcel_nbr)]*shape(parcel_nbr)] #creating a list of list. Same as cell(size(parcel_nbr)) except not empty, filled with zero
        add_vox = []
        
        for i in range(len(parcels)):
            add_mask[i]=[1 if x <= thresh else 0 for x in parcel_nbr_diff[i]] #this list comprehension gives 1 if x <= thresh, else it gives 0
            if sum(add_mask[i]) > 0:  #note summing over list not array here
                add_vox.append(parcel_nbr[i]) #list append takes only one argument
                add_vox.append(add_mask[i])
                add_vox.sort()
                array_add_vox= np.asarray(add_vox) #functions unique and diff operate on arrays, not lists; hence, converted to arrays.
                diff_list= (np.diff(array_add_vox)).tolist() 
                zero_loc = [i for i, j in enumerate(diff_list) if j == 0] #finds all positions in diff_list that equal 0.
                conflicts=np.unique(array_add_vox[zero_loc])
        
        if conflicts:
            for c in conflicts:
                #not sure of following line...
                #dists = cellfun(@(x, y) min([x(y == c) inf]), parcel_nbr_diff, parcel_nbr) 
                best_match = np.argmin(dists)
                
                for i in (np.setdiff1d(np.nonzero(np.isfinite(dists)), np.array([best_match])).tolist(): #best_match changed to array as setdiff operates on 2 arrays. Finally changed to list as it is iterable object
                    add_mask[i](parcel_nbr{i} == c) = 0 #unclear is what is meant by this statement. By "(parcel_nbr{i} == c)", is that referring to 2nd index of add_mask?

        for i in range(len(parcels)):
            #to be done for this loop.  Unclear of what is meant in expresssions like "parcel_nbr{i}(~rem_vox)"
        
    dissim = np.zeros((len(parcels),len(parcels))  
    infn = float("inf")
    
    for m in range(len(parcels)):
        for n in range(len(parcels)):
            if m == n:
                dissim[m,n] = infn
            else:
                #the following lines create same as (1:N ~= seeds(m)) & (1:N ~= seeds(n))
                ij_vector = [1]*N   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
                ij_vector[seeds[m]]=0
                ij_vector[seeds[n]]=0  

                jk=(np.arange(N)*ij_vector).tolist()  #this does something like [1,2,3,4,5] * [1,0,1,0,1], and then makes list (list because list is iterable object unlike a numpy array)
                select_vector=[x for x in jk if x!=0]  #this list comprehension only takes the values that are not equal to zero, ie [1,0,3,0,5] becomes [1,3,5]             
                dissim[m,n] = np.linalg.norm(parcel_feat[m,select_vector]-parcel_feat[n,select_vector],2))
                
    z_init = np.zeros(N)
    for i in range(len(parcels)):  
        z_init[parcels[i]] = 1
        
    parcels_adj = infn*np.ones((len(parcels),len(parcels)))
    
    for i in range(len(parcels)):
        parcel_adj[i,np.unique(z_init[adj_list[parcels[i]]])] = 1
    
    parcel_vox=cellfun(@length, parcels) #unclear...what is function being applied here...?
    
    while len(parcels) > n_clust:
        min_ind = np.argmin(dissim*parcel_adj)
        [m, n] = ind2sub(size(dissim), min_ind) #ind2sub to be implemented...
        #the following lines create other = (1:length(parcels) ~= m) & (1:length(parcels) ~= n)
        other = [1]*len(parcels)   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
        other[m]=0 
        other[n]=0  
        new_dis = (parcel_vox[m] + parcel_vox[other]) / (parcel_vox[m] + parcel_vox[n] + parcel_vox[other]) * dissim[m,other] + (parcel_vox[n] + parcel_vox[other]) / (parcel_vox[m] + parcel_vox[n] + parcel_vox[other]) * dissim[n,other] - parcel_vox[other] / (parcel_vox[m] + parcel_vox[n] + parcel_vox[other]) * dissim[m,n]
        dissim[m,other] = new_dis
        dissim[other,m] = new_dis
        
        #the following lines create (1:length(parcels) ~= n)
        jk_vector = [1]*len(parcels)   #created list as opposed to using ones((1,N)) vector because multiplication with np.arange gives a dimension for list that can't be iterated
        jk_vector[n] = 0 
        dissim = dissim(jk_vector,jk_vector)
        
        parcel_adj[m,parcel_adj[n,:]==1] = 1
        parcel_adj = parcel_adj[jk_vector,jk_vector] 
        
        parcels[m].append(parcels[n])
        parcels = parcels[jk_vector]
        
    z=np.zeros(N)
    
    for i in range(n_clust):
        z[parcels[i]] = i
    if gt_z:
        stats.NMI = CalcNMI(gt_z, z) 
    
    return(z,stats)

def CalcNMI(gt_z, z):
    N = len(gt_z)
    MI = 0;
    gt_z = gt_z.transpose() #gt_z = gt_z(:)'
    z = z.transpose() #z = z(:)'
    
    gt_p = np.zeros(np.amax(gt_z))
    H_gt = 0
    
    for i in (np.unique(gt_z)).tolist():
        gt_p[i] = np.sum(gt_z == i) / N 
        H_gt = H_gt - gt_p[i] * np.log(gt_p[i])
        
    p = np.zeros(np.amax(z))
    H = 0
    
    for j in (np.unique(z)).tolist():
        p[j]= np.sum(z == j)/N 
        H = H - p[j] * np.log(p[j])
        
    for i in (np.unique(gt_z)).tolist():
        for j in (np.unique(z)).tolist():
            joint_p = np.sum(gt_z == i and z == j) / N 
            if joint_p > 0:
                MI = MI + joint_p * np.log(joint_p / (gt_p[i]*p[j]))
                
    if MI == 0:
        NMI = 0
    else:
        NMI = MI/(np.sqrt(H*H_gt))
        
    return NMI