import numpy as np
import scipy as sp
import math as mt

def LinkageConstrained(D, adj_list):
    #what is data type of D? an array?
    if CheckSymApprox(D):
        X = D
    else:
        # D' is D.conj().transpose()
        X = np.concatenate((D,D.transpose()),axis=1)

    #dot(X,X,2) is equivalent to element squaring. tile is same as repmat    
    Qx = np.tile(np.array([(X[x])**2 for x in range(len(X))]),1,X.shape[0])
    Y = Qx + Qx.transpose()-2*np.dot(X, X.transpose())  #dot is matrix multiplication

    for i in range(len((Y.shape[0])**2)):
        Y[i,i] = 0  #sets diagonal of matrix Y to zeros...I think thats what you were trying to do..
    Y = sp.spatial.distance.squareform(Y)
    
    N = len(adj_list)
    
    A = np.zeros([N,N]) #A = false(length(adj_list)); creates a matrix of zeros or falses that is NxN in dimension
    for i in range(N):
        A[i,adj_list[i]] = 1
    connected = sp.spatial.distance.squareform(A)
    
    valid_clusts = np.ones(N) #creates array of N 1's, ie. true(N,1);
    col_limits = np.cumsum(np.arange(N-1, 0, -1))
    
    m = mt.ceil(mt.sqrt(2*Y.shape[1]))
    
    #this is probably wrong...did not understand what is meant by type "single"?
    if isinstance(Y,float32): #python types are int, long, float, complex...so not sure wha
        Z = np.zeros([m-1,3],dtype = np.float32) #float32 is same as single precision
    else:
        Z = np.zeros([m-1,3])
    
    C = np.zeros(2*m-1)
    C[0:m] = 1  #sets the first m entries as one. Selects the 0th row.
    R = np.arange(m)

    all_inds = np.arange(Y.shape[1])
    
    #all_inds is a 1D array, but connected is a matrix...???
    conn_inds = all_inds[connected] 
    
    for s in range(m-1):
        if not conn_inds:  # same as isempty(conn_inds)
            connected = np.zeros(len(connected))
            conn_inds = []
            valid_clust_inds = np.nonzero(valid_clusts) #same as find
            
            for i in valid_clust_inds: #valid_clust_inds is valid iterable
                U = valid_clusts
                U[i] = 0
                new_conns = PdistInds(i, N, U)
                connected[new_conns] = 1   #probably mistake here. what type is new_conns?
                conn_inds = np.concatenate((conn_inds, new_conns),axis=1)
            
            conn_inds = np.unique(conn_inds)  #unique is same in python

        v = np.amin(Y[conn_inds], axis=0) #find the min in each column of matrix Y.  
        k = np.argmin(Y[conn_inds], axis=0)        
    
        k = conn_inds[k.tolist()]  #unclear what is being done here..selects rows of conn_inds
        j = np.where(k <= col_limits)[0][0] #finds the first occurrance.  np.where returns a tuple, so [0][0] is first occurrence
        i = N - (col_limits[j] - k)

        Z[s,:] = np.concatenate((R[i], R[j],v),axis=1)
    
        U = valid_clusts
        i_and_j = np.concatenate((i, j),axis=1) #[i j]
        i_and_j = i_and_j.tolist()        
        U[i_and_j] = 0
        #what data type is I and J?
        I = PdistInds(i, N, U)
        J = PdistInds(j, N, U)
        
        Y[I] = ((C[R[U]]+C[R[i]])*Y[I] + (C[R[U]]+C[R[j]])*Y[J] - C[R[U]]*v)/(C[R[i]]+C[R[j]]+C[R[U]])  #element-wise computation. No dots needed in python.
        
        # update connected
        #what data type is new_conns?
        #this is same as connected(J) & ~connected(I). Can only apply ~ operator on type bool; multiplying by 1 converts to int type
        new_conns = 1*((connected[J] * 1*(~(connected[I].astype(bool)))).astype(bool))  
        connected[I] = 1*((connected[I] + new_conns).astype(bool)); #connected(I) | new_conns
        concatmat =np.concatenate((conn_inds, I(new_conns)),axis=1) #concatenate
        conn_inds = np.sort(concatmat.transpose()).transpose()  #python and matlab sort along different axis. transposing concatenated matrix, sorting, and then tranposing again appears to work
        
        U[i]=1
        J = PdistInds(j, N, U)
        
        ismem = ismembc(conn_inds,J)  #see python implementation below
        for x in ismem:    
            conn_inds[x] = []  #same as conn_inds(ismembc(conn_inds,J)) = []
       connected[J] = np.zeros(len(J))
       
       valid_clusts[j] = 0 
   
       # update m, N, R
       C[m+s] = C[R[i]] + C[R[j]]
       R[i] = m+s
       
    Z[:,2] = np.sqrt(Z[:,2]) #Z(:,3) = sqrt(Z(:,3))
    return Z

def PdistInds(row, N, valid_flags):
    if row > 0:
        N_vector=np.arange((N-2),(N-row),-1)  #Note that lowest value here is N-row and not N-row+1 as in matlab
        inds1 =np.concatenate([np.array([row-1]), (row-1)+np.cumsum(N_vector)],axis=1)  #scalar value (row-1) first has to be converted to array before concatenation 
        end = len(inds1) 
        O_vector = np.arange((inds1[end]+N-row+1),(inds1(end)+2*N-2*row)+1)
        temp_vector = np.concatenate([inds1,np.array([0])],axis=1)
        I = np.concatenate([temp_vector,O_vector],axis=1) 
    else:
        I = np.arange(N) #creates vector 0 to N-1
        
    #datatype of valid flags?    
    I = I[valid_flags];
    return I
    
    
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

def ismembc(a,b):  #python implementation of ismembc of matlab
    b_dict = {}
    for i, element in enumerate(b):
        if element not in b_dict:
            b_dict[element] = i
    return [b_dict.get(item, None) for item in a] #lookup up in hashtable is O(1)
