import numpy as np
import collections

def ddCRP(D, adj_list, init_c, gt_z, num_passes, alpha, kappa, nu, sigsq, stats_interval, verbose):
    map_z =  np.zeros(np.shape(D)[0])
    StatStruct = collections.namedtuple('Stats',['times','lp','NMI','K','z','c'])
    stats = StatStruct
    return (map_z,stats)

if __name__ == "__main__":
    res = ddCRP(np.zeros((10,10)),0,0,0,0,0,0,0,0,0,0)
    print(res[0])
