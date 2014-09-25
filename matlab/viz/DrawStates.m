function DrawStates(z)
load('../data/ACS/flowMerged.mat');
load('../data/ACS/stateData');
figure('Color',[1 1 1],'Position',[1146         646         973         587]); hold on;
for i = 1:49
    for j = 1:length(state_borders{i})
        fill(state_borders{i}{j}(:,1),state_borders{i}{j}(:,2),z(i));
    end
end
axis off;

end