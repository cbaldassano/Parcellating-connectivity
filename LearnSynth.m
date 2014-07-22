function [DC DC_K] = LearnSynth(type)

alpha=10;
kappa=0.0001;
nu=1;
sigsq = 1;
pass_limit = 20;
seeds = 12;

synth_sigsq = 0:8;
num_noise = length(synth_sigsq);
LS = zeros(num_noise,1);
NC = zeros(num_noise,1);
WC = zeros(num_noise,1);
RG = zeros(num_noise,1);
DC = zeros(num_noise,1);
DC_stats = cell(num_noise, seeds);
DC_K = zeros(num_noise,1);

if (strcmp(type, 'square'))
    n_clust = 9;
elseif (strcmp(type, 'stripes'))
    n_clust = 6;
elseif (strcmp(type, 'face'))
    n_clust = 5;
end

subject = 'synth';
for noise_ind = 1:num_noise
    disp(['Noise level ' num2str(noise_ind)]);
    experiment = [type '_' num2str(synth_sigsq(noise_ind))];
    
    loaded = load(['../data/' subject '/' experiment '.mat']);
    D = loaded.D;
    adj_list = loaded.adj_list;
    coords = loaded.coords;
    gt_z = loaded.z;
    clear loaded;
    
%     disp('   Local Similarity');
%     z = LocalSimilarity(D, adj_list, n_clust);
%     LS(noise_ind) = CalcNMI(gt_z, z);
%     
%     disp('   Ward Clustering');
%     z = WardClustering(D, adj_list, n_clust);
%     WC(noise_ind) = CalcNMI(gt_z, z);
%     
%     disp('   Normalized Cut');
%     z = NormCut(D, adj_list, n_clust);
%     NC(noise_ind) = CalcNMI(gt_z, z);
    
%     disp('   Region Growing');
%     addpath('reggrow');
%     addpath('reggrow/hlpfunc');
%     max_adj = max(cellfun(@length,adj_list));
%     padded_adj_list = cellfun(@(x) [x zeros(1,max_adj-length(x))], ...
%         adj_list, 'UniformOutput', false);
%     z = t_reggrow_main(D, cell2mat(padded_adj_list), 3*[coords zeros(size(D,1),1)], n_clust);
%     RG(noise_ind) = CalcNMI(gt_z, z);
     
    disp('   ddCRP');
    D = NormalizeConn(D);
    [~, Z] = WardClustering(D, adj_list, 1:20);
    parfor seed = 1:seeds
        %disp(['     Seed ' num2str(seed)]);
        rng(seed);
        [~,stats] = InitializeAndRunddCRP(Z, D, adj_list, 1:20, alpha, kappa, nu, sigsq, pass_limit, gt_z, 0);
        DC_stats{noise_ind,seed} = stats;
    end
    
    max_lp = -inf;
    for seed = 1:seeds
        if (DC_stats{noise_ind, seed}.lp(end) > max_lp)
            max_lp = DC_stats{noise_ind, seed}.lp(end);
            max_lp_seed = seed;
        end
    end
    DC(noise_ind) = DC_stats{noise_ind, max_lp_seed}.NMI(end);
    DC_K(noise_ind) = DC_stats{noise_ind, max_lp_seed}.K(end);
end

%save(['../output/synth/' type '.mat'], 'synth_sigsq','LS', 'NC', 'WC', 'RG', 'DC','DC_stats');
end

