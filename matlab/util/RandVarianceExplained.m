function rand_var = RandVarianceExplained(K, num_seeds)

load('../data/Q1/100307/PPA.mat');
rand_var = zeros(num_seeds,1);
for i =1:num_seeds
    init_z = labels;
    init_z(labels>=9) = 0;
    rand_z = GenerateRandomClustering(adj_list, K, init_z');
    rand_var(i) = CalcVarianceExplained(D, rand_z');
end

end

