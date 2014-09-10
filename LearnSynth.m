function [LS NC WC RG DC DC_K RC] = LearnSynth(type)

rng(1);

alpha=10;
kappa=0.0001;
nu=1;
sigsq = 0.01;
pass_limit = 30;
repeats = 20;

synth_sig = linspace(0,9,10);
num_noise = length(synth_sig);
LS = zeros(num_noise,repeats);
NC = zeros(num_noise,repeats);
WC = zeros(num_noise,repeats);
RG = zeros(num_noise,repeats);
DC = zeros(num_noise,repeats);
DC_K = zeros(num_noise,repeats);
RC = zeros(num_noise, repeats);

for rep = 1:repeats
    for noise_ind = 1:num_noise
        disp(['Noise level ' num2str(noise_ind)]);
        [D adj_list gt_z coords] = GenerateSynthData(type, synth_sig(noise_ind));

        disp('   ddCRP');
        D = NormalizeConn(D);
        [~, Z] = WardClustering(D, adj_list, 1);
        [~, stats] = InitializeAndRunddCRP(Z, D, adj_list, 1:20, alpha, kappa, nu, sigsq, pass_limit, [], gt_z, 0);
        DC(noise_ind, rep) = stats.NMI(end);
        DC_K(noise_ind, rep) = stats.K(end);
        
        n_clust = DC_K(noise_ind, rep);
        
        disp('   Local Similarity');
        z = LocalSimilarity(D, adj_list, n_clust);
        LS(noise_ind, rep) = CalcNMI(gt_z, z);
 
        disp('   Ward Clustering');
        z = WardClustering(D, adj_list, n_clust);
        WC(noise_ind, rep) = CalcNMI(gt_z, z);

        disp('   Normalized Cut');
        z = NormCut(D, adj_list, n_clust);
        NC(noise_ind, rep) = CalcNMI(gt_z, z);

        disp('   Region Growing');
        addpath('reggrow');
        addpath('reggrow/hlpfunc');
        max_adj = max(cellfun(@length,adj_list));
        padded_adj_list = cellfun(@(x) [x zeros(1,max_adj-length(x))], ...
            adj_list, 'UniformOutput', false);
        z = t_reggrow_main(D, cell2mat(padded_adj_list), 3*[coords zeros(size(D,1),1)], n_clust);
        RG(noise_ind, rep) = CalcNMI(gt_z, z);

        disp('   Random Clustering');
        z = GenerateRandomClustering(adj_list, n_clust);
        RC(noise_ind, rep) = CalcNMI(gt_z, z);
    end
end
end

