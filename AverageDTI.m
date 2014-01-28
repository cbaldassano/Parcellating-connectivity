function AverageDTI()

D = zeros(59412);
D = single(D);
for subj = {'100408', '105216', '106319', '111514', '112819', '101915', '102816', '106016', '111009', '111716'}
    loaded = load(['../data/Q3/' subj{1} '/full.mat']);
    D = D + loaded.D/10;
end
save('../data/Q3/group/full.mat','D', '-v7.3');
end

