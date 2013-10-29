function LearnSeeds(sigsq)

num_seeds = 5;
map_z = cell(num_seeds,1);
stats = cell(num_seeds,1);
parfor seed = 1:num_seeds
    rng(seed); 
    [map_z{seed} stats{seed}] = ddCRP('unrelated40','full',5,10,0.0001,1,sigsq,1000,1);
end

save(['../output/unrelated40/full_1pass_' num2str(sigsq) '.mat'],'map_z','stats');


end

