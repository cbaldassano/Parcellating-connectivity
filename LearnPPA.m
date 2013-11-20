function LearnPPA()
subj = {'100307' '103515' '103818' '111312' '117122' '118932' '119833' '120212' '125525' ...
        '128632' '130013' '137128' '138231' '142828' '143325' '144226' '149337' '156637' '159239' '161731'};
sigsq = 1:10;
%num_seeds = 5;
%map_z = cell(length(subj), length(sigsq), num_seeds);
%stats = cell(length(subj), length(sigsq), num_seeds);
map_z = cell(length(subj), length(sigsq));
stats = cell(length(subj), length(sigsq));
for s = 1:length(subj)
    disp(subj{s});
    parfor ss = 1:length(sigsq)
        %parfor seed = 1:num_seeds
        %    rng(seed); 
        %    [map_z{s,ss,seed} stats{s,ss,seed}] = ddCRP(['Q1/' subj{s}],'PPA',20,10,0.0001,1,sigsq(ss),100,1);
        %end
        [map_z{s,ss} stats{s,ss}] = ddCRP(['Q1/' subj{s}],'PPA',10,10,0.0001,1,sigsq(ss),100,1);
    end
end

save('../output/PPA/allsubj.mat','map_z','stats');


end

