function [bold, D] = LoadSubjRFMRI(datapath)

disp(['Loading from ' datapath '...']);
metric_files = {'lr1', 'lr2', 'rl1', 'rl2'};
hems = {'left', 'right'};

bold = cell(2,1);
bold{1} = zeros(32492, 1200*4, 'single');
bold{2} = zeros(32492, 1200*4, 'single');


for hem = 1:2
    for i = 1:length(metric_files)
        fid = fopen(fullfile(datapath,[metric_files{i} '_' hems{hem} '.metric']));
        xmlstr = fread(fid,'*char')';
        fclose(fid);
        block_starts = strfind(xmlstr, '<Data>');
        block_starts = block_starts + 6;
        block_ends = strfind(xmlstr, '</Data>');
        block_ends = block_ends - 1;
        
        for t = 1:1200
            bold{hem}(:,1200*(i-1) + t) = typecast(dunzip(base64decode( ...
                xmlstr(block_starts(t):block_ends(t)))), 'single');
        end
        bold{hem}(:,(1200*(i-1)+1):(1200*i)) = bsxfun(@minus, bold{hem}(:,(1200*(i-1)+1):(1200*i)), mean(bold{hem}(:,(1200*(i-1)+1):(1200*i)),2));
    end
    
    % Remove invalid voxels
    valid_vox = any(bold{hem},2);
    bold{hem} = bold{hem}(valid_vox,:);
end

bold = cell2mat(bold);
if (nargout > 1)
    D = atanh(corr(bold'));
end

end