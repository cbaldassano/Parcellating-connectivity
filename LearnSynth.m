function LearnSynth(type)

alpha=10;
kappa=0.0001;
nu=1;
sigsq = 1;
pass_limit = 10;
seeds = 10;

synth_sigsq = 0:8;
num_noise = length(synth_sigsq);
LS = zeros(num_noise,1);
NC = zeros(num_noise,1);
RG = zeros(num_noise,1);
DC = zeros(num_noise,1);
DC_stats = cell(num_noise, seeds);

if (strcmp(type, 'square'))
    n_clust = 9;
elseif (strcmp(type, 'stripes'))
    n_clust = 6;
elseif (strcmp(type, 'face'))
    n_clust = 5;
end

for noise_ind = 1:num_noise
    disp(['Noise level ' num2str(noise_ind)]);
    
    disp('   Local Similarity');
    [~, stats] = ...
        LocalSimilarity('synth', ...
        [type '_' num2str(synth_sigsq(noise_ind))], n_clust);
    LS(noise_ind) = stats.NMI;
    
    disp('   Normalized Cut');
    [~, stats] = ...
        NormCut('synth', ...
        [type '_' num2str(synth_sigsq(noise_ind))], n_clust);
    NC(noise_ind) = stats.NMI;
    
    disp('   Region Growing');
    [~, stats] = ...
        RegionGrowing('synth', ...
        [type '_' num2str(synth_sigsq(noise_ind))], n_clust);
    RG(noise_ind) = stats.NMI;
    
    disp('   ddCRP');
    parfor seed = 1:seeds
        disp(['     Seed ' num2str(seed)]);
        rng(seed);
        [~,stats] = ddCRP('synth', ...
            [type '_' num2str(synth_sigsq(noise_ind))], pass_limit, ...
            alpha, kappa, nu, sigsq, 1000, 0);
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
end

save(['../output/synth/' type '.mat'], 'synth_sigsq','LS', 'NC', 'RG', 'DC','DC_stats');
end

