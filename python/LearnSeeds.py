import numpy as np
import ddCRP

def LearnSeeds(sigsq):
    num_seeds = 5
    map_z = np.empty(num_seeds,dtype=object)
    stats = np.empty(num_seeds,dtype=object)

    for seed in range(num_seeds):   #matlab implementation uses parfor...
        np.random.seed(seed)
        map_z[seed], stats[seed] = ddCRP('unrelated40','full',10,10,0.0001,1,sigsq,1000,1);
    
    #save not implemented...