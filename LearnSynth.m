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
DC = cell(num_noise, seeds);

if (strcmp(type, 'square'))
    n_clust = 9;
elseif (strcmp(type, 'stripes'))
    n_clust = 6;
elseif (strcmp(type, 'face'))
    n_clust = 5;
end

for noise_ind = 1:num_noise;
    [~, stats] = ...
        LocalSimilarity('synth', ...
        [type '_' num2str(synth_sigsq(noise_ind))], n_clust);
    LS(noise_ind) = stats.NMI;
    
    [~, stats] = ...
        NormCut('synth', ...
        [type '_' num2str(synth_sigsq(noise_ind))], n_clust);
    NC(noise_ind) = stats.NMI;
    
    [~, stats] = ...
        RegionGrowing('synth', ...
        [type '_' num2str(synth_sigsq(noise_ind))], n_clust);
    RG(noise_ind) = stats.NMI;
    
%     for seed = 1:seeds
%         rng(seed);
%         [~,stats] = ddCRP('synth', ...
%             [type '_' num2str(synth_sigsq(noise_ind))], pass_limit, ...
%             alpha, kappa, nu, sigsq);
%         DC{noise_ind,seed} = stats.NMI;
%     end
end

save(['../output/synth/' type '.mat'], 'LS', 'NC', 'RG', 'DC');
end

