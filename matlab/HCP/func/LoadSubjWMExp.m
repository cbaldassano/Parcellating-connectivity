function D = LoadSubjWMExp(datapath)

files = {'tfMRI_WM_LR/LR_left.metric', 'tfMRI_WM_LR/LR_right.metric';'tfMRI_WM_RL/RL_left.metric', 'tfMRI_WM_RL/RL_right.metric'};

TRs = 405;

bold = cell(2,1);
bold{1} = zeros(32492, TRs*2);
bold{2} = zeros(32492, TRs*2);


for hem = 1:2
    for i = 1:size(files,1)
        fid = fopen(fullfile(datapath,files{i,hem}));
        xmlstr = fread(fid,'*char')';
        fclose(fid);
        block_starts = strfind(xmlstr, '<Data>');
        block_starts = block_starts + 6;
        block_ends = strfind(xmlstr, '</Data>');
        block_ends = block_ends - 1;
        
        for t = 1:TRs
            bold{hem}(:,TRs*(i-1) + t) = typecast(dunzip(base64decode(xmlstr(block_starts(t):block_ends(t)))), 'single');
        end
        bold{hem}(:,(TRs*(i-1)+1):(TRs*i)) = bsxfun(@minus, bold{hem}(:,(TRs*(i-1)+1):(TRs*i)), mean(bold{hem}(:,(TRs*(i-1)+1):(TRs*i)),2));
    end
    
    % Remove invalid voxels
    valid_vox = any(bold{hem},2);
    bold{hem} = bold{hem}(valid_vox,:);
end

D = cell2mat(bold);

end