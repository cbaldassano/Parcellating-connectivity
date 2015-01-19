function DrawCounties(z)
load('../../../data/ACS/flowMerged.mat');
load('../../../data/ACS/stateData');
figure('Color',[1 1 1],'Position',[1146         646         973         587]); hold on;
for i = 1:2594
    for j = 1:length(borders{i})
        fill(borders{i}{j}(:,1),borders{i}{j}(:,2),z(i),'EdgeColor','none');
    end
end
for i = 1:49
    for j = 1:length(state_borders{i})
        plot(state_borders{i}{j}(:,1),state_borders{i}{j}(:,2),'k');
    end
end
axis off;

K = max(z);
z_adj = cell(K,1);
for i = 1:K
    z_adj{i} = unique(z(setdiff([adj_list{z == i}],find(z == i))));
end
palette = PTPalette(12);
colors = repmat(palette,ceil(K/12),1);
colors = colors(1:K,:);
no_conflict = false;
while ~no_conflict
    no_conflict = true;
    for i = 11:K
        if (ismember(colors(i,:),colors(z_adj{i},:),'rows'));
            colors(i,:) = palette(randi(12),:);
            no_conflict = false;
        end
    end
end
colormap(colors);
end