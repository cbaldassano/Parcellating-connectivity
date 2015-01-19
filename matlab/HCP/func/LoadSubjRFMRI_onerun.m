function [bold, D] = LoadSubjRFMRI_onerun(datapath, metric)

hems = {'left', 'right'};

bold = cell(2,1);
bold{1} = zeros(32492, 1200, 'single');
bold{2} = zeros(32492, 1200, 'single');


for hem = 1:2
    fid = fopen(fullfile(datapath,[metric '_' hems{hem} '.metric']));
    xmlstr = fread(fid,'*char')';
    fclose(fid);
    block_starts = strfind(xmlstr, '<Data>');
    block_starts = block_starts + 6;
    block_ends = strfind(xmlstr, '</Data>');
    block_ends = block_ends - 1;

    for t = 1:1200
        bold{hem}(:,t) = typecast(dunzip(base64decode( ...
            xmlstr(block_starts(t):block_ends(t)))), 'single');
    end
    bold{hem} = bsxfun(@minus, bold{hem}, mean(bold{hem},2));
    
    % Remove invalid voxels
    valid_vox = any(bold{hem},2);
    bold{hem} = bold{hem}(valid_vox,:);
end

bold = cell2mat(bold);
if (nargout > 1)
    D = atanh(corr(bold'));
end

end