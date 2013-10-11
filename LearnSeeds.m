function LearnSeeds(sigsq)

map_z = cell(10,1);
stats = cell(10,1);
parfor seed = 1:10
    rng(seed); 
    [map_z{seed} stats{seed}] = ddCRP('unrelated40','small',20,10,0.0001,1,sigsq,0);
end

save(['../output/unrelated40/small_' num2str(sigsq) '.mat'],'map_z','stats');


end

