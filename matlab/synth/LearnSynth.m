% Computes a parcellation of synthetic data at different noise
%   levels, using Ward Clustering and our method based on the ddCRP. Each
%   parcellation is evaluated based on its Normalized Mututal Information
%   with the ground truth. The input "type"={'square','stripes','face'}
%   determines the underlying ground truth parcellation.
function [LS NC WC RG DC DC_K RC] = LearnSynth(type,wardsize)

rng(1); % For repeatability

% Hyperparameters
alpha=10;
kappa=0.0001;
nu=1;
sigsq = 0.01;
pass_limit = 100;%30;
repeats = 2;   % Number of times to repeat experiments

synth_sig = linspace(0,9,10);   % Noise levels to try
num_noise = length(synth_sig);
LS = zeros(num_noise,repeats);
NC = zeros(num_noise,repeats);
WC = zeros(num_noise,repeats);
RG = zeros(num_noise,repeats);
DC = zeros(num_noise,repeats);
DC_K = zeros(num_noise,repeats);
RC = zeros(num_noise, repeats);

for rep = 1:repeats
    disp(['Repeat # ' num2str(rep)]);
    parfor noise_ind = 1:num_noise
        [D adj_list gt_z coords] = GenerateShapeData(synth_sig(noise_ind)); %GenerateSynthData(type, synth_sig(noise_ind));

        D = NormalizeConn(D);
        [~, Z] = WardClustering(D, adj_list, 1);
        [~, stats] = InitializeAndRunddCRP(Z, D, adj_list, wardsize, alpha, kappa, nu, sigsq, pass_limit, [], gt_z, 0); %InitializeAndRunddCRP(Z, D, adj_list, 1:20, alpha, kappa, nu, sigsq, pass_limit, [], gt_z, 0);
        DC(noise_ind, rep) = stats.NMI(end);
        DC_K(noise_ind, rep) = stats.K(end);
        
%         n_clust = DC_K(noise_ind, rep);
%         
%         z = LocalSimilarity(D, adj_list, n_clust);
%         LS(noise_ind, rep) = CalcNMI(gt_z, z);
%  
%         z = WardClustering(D, adj_list, n_clust);
%         WC(noise_ind, rep) = CalcNMI(gt_z, z);
% 
%         z = NormCut(D, adj_list, n_clust);
%         NC(noise_ind, rep) = CalcNMI(gt_z, z);
% 
%         max_adj = max(cellfun(@length,adj_list));
%         padded_adj_list = cellfun(@(x) [x zeros(1,max_adj-length(x))], ...
%             adj_list, 'UniformOutput', false);
%         z = t_reggrow_main(D, cell2mat(padded_adj_list), 3*[coords zeros(size(D,1),1)], n_clust);
%         RG(noise_ind, rep) = CalcNMI(gt_z, z);
% 
%         z = GenerateRandomClustering(adj_list, n_clust);
%         RC(noise_ind, rep) = CalcNMI(gt_z, z);
    end
end
end

