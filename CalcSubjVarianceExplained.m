function [subj_var] = CalcSubjVarianceExplained(z)
subj = {'120515','113922','111413','109325','109123','108828','108525','108323','108121','107422','107321','106521','105014','104820','103111','102311','102008','101309','101107','101006'};

subj_var = zeros(length(subj),1);
for i = randperm(length(subj))
    load(['../data/S500/' subj{i} '/bold']);
    bold = cast(bold, 'double');
    bold = bsxfun(@minus, bold, mean(bold,2));
    bold = bsxfun(@rdivide, bold, sqrt(sum(bold.*bold,2)));
    N = size(bold,1);
    
    % Calc base var
    if (exist(['../data/S500/' subj{i} '/basevar.mat'],'file'))
        load(['../data/S500/' subj{i} '/basevar.mat']);
    else
        base_sum = 0;
        for v1 = 1:N
            %if (mod(v1,100)==1)
            %    disp([num2str(v1/N*100) '%']);
            %end
            corrs = bold*(bold(v1,:)');
            base_sum = base_sum + sum(atanh(corrs(1:N > v1)));
        end
        base_mean = base_sum/(N*(N-1)/2);
        base_var = 0;
        for v1 = 1:N
            %if (mod(v1,100)==1)
            %    disp([num2str(v1/N*100) '%']);
            %end
            corrs = bold*(bold(v1,:)');
            base_var = base_var + sum((base_mean - atanh(corrs(1:N > v1))).^2);
        end
        save(['../data/S500/' subj{i} '/basevar.mat'],'base_var');
        disp(['Saved ' subj{i} ' base_var']);
    end
    
    z = z(:)';
    [sorted_z, sorted_i] = sort(z);
    bins = mat2cell(sorted_i, 1, diff(find(diff([0 sorted_z (max(z)+1)]))));
    
    sum_var = 0;
    for v1=1:length(bins)
        %disp(['Bins: ' num2str(v1/length(bins)*100) '%']);
        for v2=v1:length(bins)
            x = atanh(bold(bins{v1},:)*(bold(bins{v2},:)'));
            if (v1==v2)
                x = x(~tril(ones(size(x))));
            else
                x = x(:);
            end
            sum_var = sum_var + sum((x - mean(x)).^2);
        end
    end
    subj_var(i) = 1 - sum_var / base_var;
    disp(['Subj ' subj{i} ': ' num2str(subj_var(i))]);
end
end
