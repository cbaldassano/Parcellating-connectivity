import numpy as np
import collections
import WardClustering
import StatsUtil
import InitializeAndRunddCRP as initdd

def LearnSynth(type):

    np.random.seed(1)
    max_noise = 10;
    repeats = 20;

    # ddCRP hyperparameters
    alpha = 10;
    kappa = 0.0001;
    nu = 1;
    sigsq = 0.01;
    pass_limit = 30;
    
    WC = np.zeros((max_noise,repeats))
    DC = np.zeros((max_noise,repeats))
    DC_K = np.zeros((max_noise,repeats))

    for rep in range(repeats):
        print('Repeat #' + str(rep))
        for noise_sig in range(max_noise):
            print('   Noise level: ' + str(noise_sig))
            synth = GenerateSynthData(type, noise_sig)
            D = NormalizeConn(synth.D)

            # ddCRP
            Z = WardClustering.ClusterTree(D, synth.adj_list)
            dd_results = initdd.InitializeAndRun(Z, D, synth.adj_list, range(1,21), alpha, kappa, nu, sigsq, pass_limit, synth.z, 0)
            DC[noise_sig,rep] = dd_results[1].NMI[-1]
            DC_K[noise_sig,rep] = dd_results[1].K[-1]

            n_clust = DC_K[noise_sig,rep]

            # Ward Clustering
            WC[noise_sig, rep] = StatsUtil.NMI(synth.z, WardClustering.Cluster(Z, n_clust))

    return (WC,DC,DC_K)

# Generate synthetic dataset (connectivity, adjacency, and coordinates) at given noise level
def GenerateSynthData(type, sig):
    sqrtN = 18
    
    SynthData = collections.namedtuple('SynthData',['D','adj_list','z','coords'])
    coords = np.zeros((sqrtN**2,2))
    adj_list = np.empty(sqrtN**2, dtype=object)
    for r in range(0, sqrtN):
        for c in range(0, sqrtN):
            currVox = c + r*sqrtN
            coords[currVox,:] = [r, c]
            curr_adj = []
            if r > 0:
                curr_adj.append(c + (r-1)*sqrtN)
            if r < (sqrtN-1):
                curr_adj.append(c + (r+1)*sqrtN)
            if c > 0:
                curr_adj.append((c-1) + r*sqrtN)
            if c < (sqrtN-1):
                curr_adj.append((c+1) + r*sqrtN)
            adj_list[currVox] = np.array(curr_adj)
    
    if type == 'square':
        z = np.array([
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8])
    elif type == 'stripes':
        z = np.array([
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        0,0,0,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,
        0,0,0,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,
        0,0,0,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,
        0,0,0,1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,
        0,0,0,1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,
        0,0,0,1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,
        0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,4,4,4,
        0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,4,4,4,
        0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,4,4,4,
        0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,
        0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,
        0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5])
    elif type == 'face':
        z = np.array([
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        0,0,0,0,0,0,3,3,3,3,3,3,6,6,6,6,6,6,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        1,1,1,1,1,1,4,4,4,4,4,4,7,7,7,7,7,7,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8,
        2,2,2,2,2,2,5,5,5,5,5,5,8,8,8,8,8,8])
    
    D = GenConnectivity(z, sig)
    
    synth = SynthData(D, adj_list, z, coords)
    return synth

# Generate synthetic connectivity matrix at given noise level
def GenConnectivity(z, sig):
    N = len(z)
    K = len(np.unique(z))
    
    A = np.random.normal(size=(K,K))
    
    D = np.zeros((N,N))
    for v1 in range(0,N):
        for v2 in range(0,N):
            if v1 != v2:
                D[v1,v2] = sig*np.random.normal() + A[z[v1],z[v2]]
    
    return D
    

# Normalize connectivity matrix to have zero mean and unit variance
def NormalizeConn(D):
    D = D.astype('float64')
    N = D.shape[0]
    off_diags = np.logical_not(np.eye(N,dtype='bool'))
    D = D - D[off_diags].mean()
    D = D/D[off_diags].std()
    np.fill_diagonal(D, 0)

    D = D.astype('float32')
    return D

if __name__ == "__main__":
    NMIs = LearnSynth('stripes');
    print('WC: ' + NMIs[0])
    print('DC: ' + NMIs[1])
