function MAP4=t_grow_regions_from_seeds(Data, OUT, OUT2, Mask, NB, USE_ONLY_CENTRE ,use_dist)

MaskedOut=find(Mask==0);
fz=find(OUT(:,1)==0);
CENTRES=OUT2;
CENTRES(fz)=0;
fX=find(CENTRES);
lfX=length(fX);
D3=zeros(size(Data,1),lfX);
ctr=1;
for ctr1=fX'
    if USE_ONLY_CENTRE
         D3(:,ctr)=Data(:,ctr1);
         ctr=ctr+1;
    else
         nbrs           = NB(ctr1);
         nbrs           = unique(nbrs(:));
         nbrs(nbrs==0)  = [];
         nbrs           = union(nbrs,ctr1);
         nbrs           = intersect(nbrs,find(Mask));
         D3(:,ctr)      = mean(Data(:,nbrs),2);
         D3(:,ctr)      = D3(:,ctr)-mean(D3(:,ctr));
         STD            = std(D3(:,ctr));
         if STD>1e-12
             D3(:,ctr)=D3(:,ctr)/STD;
         else
             error('ERROR:STD of timecourse is too small, make sure this does not happen.')
         end
    end
end
% CB computing correlation explicitly
% M2=D3'*Data./(size(Data,1)-1);
% distMAT=0.5+0.5*M2;
distMAT = 0.5 + 0.5*corr(D3,Data);

ldat=size(Data,2);
CLUSTS=spalloc(lfX,ldat,ldat);
cl=ones(lfX,1);
for ctr=1:lfX
    CLUSTS(ctr,1)=fX(ctr);
end
uC=unique(fX);

% INITIALISE
CMAT=sparse(lfX,ldat); 
MAX_N=zeros(lfX,1);
MAX_C=zeros(lfX,1);

%lfX
for ctr2=1:lfX
    nbrs=NB(CLUSTS(ctr2,1:cl(ctr2)),:);
    nbrs=nbrs(:)';
    nbrs(nbrs==0)=[];
    nbrs=setdiff(nbrs,union(MaskedOut,[CLUSTS(:)]));
    CMAT(ctr2,nbrs)=distMAT(ctr2,nbrs);
end

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SLOW CENTRE BASED               
if sum(strcmp(use_dist,{'centre';'mean';'max';'min';'ward'}))
    [a b]=max(CMAT,[],1);
    sCMAT=sparse(b, 1:length(b),a, lfX,size(Data,2));
    [MAX_C,MAX_N]=max(sCMAT,[],2);
    RUN=1;
    ctr=1;
    CUT=0.01;
    while RUN
        if ctr/size(Data,2)>CUT
            CUT=CUT+0.01;
            ctr/size(Data,2);
        end
        [val idx] = max(MAX_C);
        if MAX_C(idx)>0
            NEW_IDX=MAX_N(idx);
            CLUSTS(idx,cl(idx)+1:cl(idx)+1)
            uC=[uC NEW_IDX];
            cl(idx)=cl(idx)+1;
            % UPDATE DISTS
            if strcmp(use_dist,'max')
                distMAT(idx,:)=(max([distMAT(idx,:)' Data'*Data(:,idx)],[],2))';
            elseif strcmp(use_dist,'min')
                distMAT(idx,:)=(min([distMAT(idx,:)' Data'*Data(:,idx)],[],2))';
            elseif strcmp(use_dist,'mean')
                distMAT(idx,:)=(mean([distMAT(idx,:)' Data'*Data(:,idx)],2))';   
            elseif strcmp(use_dist,'ward')
                n=cl(idx)-1;
                distMAT(idx,:)=(1/(n+1)*sum([n*distMAT(idx,:)' Data'*Data(:,idx)],2))'; 
            elseif strcmp(use_dist,'centre')
                %NOTHING HAPPENS
            end

            % update CMAT
            % First set column for selected index to zero
            CMAT(:,NEW_IDX)=0;
            % Add elements for new neighbours of just grown region
            %fCMAT=find(CMAT(idx,:));
            nbrs=NB(CLUSTS(idx,1:cl(idx)),:);
            nbrs=nbrs(:)';
            nbrs(nbrs==0)=[];
            nbrs=setdiff(nbrs,union(MaskedOut,uC));
            CMAT(idx,nbrs)=distMAT(idx,nbrs);
        else
            RUN=0;
        end
        [a b]=max(CMAT,[],1);
        sCMAT=sparse(b, 1:length(b),a, lfX,size(Data,2));
        [MAX_C,MAX_N]=max(sCMAT,[],2);
        ctr=ctr+1;
        %length(unique(MAX_N(MAX_N~=0)))/length(MAX_N(MAX_N~=0))
    end
    %
    MAP4=zeros(size(Data,2),1);
    for ctr=1:size(CLUSTS,1);
        MAP4(CLUSTS(ctr,1:cl(ctr)))=ctr;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FAST CENTRE BASED    
elseif strcmp(use_dist,'fast')
    val=max(max(CMAT));
    RUN=1;
    ctr=0;
    CUT=0.01;
    time_run=0;
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    % relaxation parameter
    sc=0.9;
    tic
    while RUN
       % if ctr/size(Data,2)>CUT
       %     CUT=CUT+0.01;

            percent_done = ctr/size(Data,2);
            if ctr>0
                time_run     = time_run+toc/60;
                approx_time_to_go   = time_run/percent_done-time_run;
            end
            tic
        %end
        Ic=[];
        if val>0
            fl      = find(CMAT>=sc*val);
            [Ir Ic] = ind2sub(size(CMAT),fl);
            % See if there are several vertices that are to be assigned to different regions, if so, pick largest and delete the others 
            [u uid] = unique(Ic);
            cid     = setdiff(1:length(Ic),uid);
            if ~isempty(cid)
                cv      = Ic(cid); 
                cv      = unique(cv);
                for ctr5=1:length(cv)
                   fc       = find(Ic==cv(ctr5));
                   cvals    = CMAT(sub2ind(size(CMAT),Ir(fc),Ic(fc)));
                   %Find largest
                   [~, ci]  = max(cvals);
                   % Delete others
                   fc(ci)   = [];
                   Ir(fc)   = [];
                   Ic(fc)   = [];
                   fl(fc)   = [];
                end
            end
            uI=unique(Ir);
            for ctr5=1:length(uI)
                idx=uI(ctr5);
                fc=find(Ir==idx);
                CLUSTS(idx,cl(idx)+1:cl(idx)+length(fc))=[Ic(fc)]';
                uC=[uC; Ic(fc)];
                cl(idx)=cl(idx)+length(fc);
                % update CMAT
                % First set column for selected index to zero
                CMAT(:,Ic(fc))=0;
                % Add elements for new neighbours of just grown region
                % fCMAT=find(CMAT(idx,:));
                nbrs=NB(CLUSTS(idx,1:cl(idx)),:);
                nbrs=nbrs(:)';
                nbrs(nbrs==0)=[];
                %nbrs=setdiff(nbrs,union(MaskedOut,union([CLUSTS(:)],fCMAT)));

                nbrs=setdiff(nbrs,[MaskedOut; uC]);

                CMAT(idx,nbrs)=distMAT(idx,nbrs);
            end
        else
            RUN=0;
        end
        val=max(max(CMAT));
        ctr=ctr+length(Ic);
    end
    MAP4=zeros(size(Data,2),1);
    for ctr=1:size(CLUSTS,1);
        MAP4(CLUSTS(ctr,1:cl(ctr)))=ctr;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FAST MEAN BASED   
elseif strcmp(use_dist,'fast_mean')

    val=max(max(CMAT));
    RUN=1;
    ctr=0;
    CUT=0.01;
    time_run=0;
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    % relaxation parameter
    sc=0.9;
    tic
    while RUN
       % if ctr/size(Data,2)>CUT
       %     CUT=CUT+0.01;

        percent_done = ctr/size(Data,2)
        if ctr>0
            time_run     = time_run+toc/60
            approx_time_to_go   = time_run/percent_done-time_run
        end
        tic
        %end
        if val>0
            fl      = find(CMAT>=sc*val);
            [Ir Ic] = ind2sub(size(CMAT),fl);
            % See if there are several vertices that are to be assigned to different regions, if so, pick largest and delete the others 
            [u uid] = unique(Ic);
            cid     = setdiff(1:length(Ic),uid);
            if ~isempty(cid)
                cv      = Ic(cid); 
                cv      = unique(cv);
                for ctr5=1:length(cv)
                   fc       = find(Ic==cv(ctr5));
                   cvals    = CMAT(sub2ind(size(CMAT),Ir(fc),Ic(fc)));
                   %Find largest
                   [~, ci]  = max(cvals);
                   % Delete others
                   fc(ci)   = [];
                   Ir(fc)   = [];
                   Ic(fc)   = [];
                   fl(fc)   = [];
                end
            end
            uI=unique(Ir);
            for ctr5=1:length(uI)
                idx=uI(ctr5);
                fc=find(Ir==idx);
                CLUSTS(idx,cl(idx)+1:cl(idx)+length(fc))=[Ic(fc)]';
                uC=[uC; Ic(fc)];
                cl(idx)=cl(idx)+length(fc);

                % update CMAT
                % First set column for selected index to zero
                CMAT(:,Ic(fc))=0;
                % Add elements for new neighbours of just grown region
                % fCMAT=find(CMAT(idx,:));
                nbrs=NB(CLUSTS(idx,1:cl(idx)),:);
                nbrs=nbrs(:)';
                nbrs(nbrs==0)=[];
                %nbrs=setdiff(nbrs,union(MaskedOut,union([CLUSTS(:)],fCMAT)));
                nbrs=setdiff(nbrs,[MaskedOut; uC]);

                inds          = CLUSTS(idx,1:cl(idx));

                D3(:,idx)     = mean(Data(:,inds),2);
                D3(:,idx)     = D3(:,idx)-mean(D3(:,idx));
                STD           = std(D3(:,idx));
                if STD>1e-12
                    D3(:,idx)=D3(:,idx)/STD;
                else
                    error('ERROR:STD of timecourse is too small, make sure this does not happen.')
                end
                m2=(D3(:,idx)'*Data(:,nbrs))./(size(Data,1)-1);
                %distMAT(idx,nbrs)=0.5+0.5*m2;
                %CMAT(idx,nbrs)=distMAT(idx,nbrs);
                CMAT(idx,:)=0;
                CMAT(idx,nbrs)=0.5+0.5*m2;
            end
        else
            RUN=0;
        end
        val=max(max(CMAT));
        ctr=ctr+length(Ic);
    end
    MAP4=zeros(size(Data,2),1);
    for ctr=1:size(CLUSTS,1);
        MAP4(CLUSTS(ctr,1:cl(ctr)))=ctr;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FAST MIN BASED   
elseif strcmp(use_dist,'fast_min')

    val=max(max(CMAT));
    RUN=1;
    ctr=0;
    CUT=0.01;
    time_run=0;
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    % relaxation parameter
    sc=0.9;
    tic
    while RUN
       % if ctr/size(Data,2)>CUT
       %     CUT=CUT+0.01;

        percent_done = ctr/size(Data,2)
        if ctr>0
            time_run     = time_run+toc/60
            approx_time_to_go   = time_run/percent_done-time_run
        end
        tic
        %end
        if val>0
            fl      = find(CMAT>=sc*val);
            [Ir Ic] = ind2sub(size(CMAT),fl);
            % See if there are several vertices that are to be assigned to different regions, if so, pick largest and delete the others 
            [u uid] = unique(Ic);
            cid     = setdiff(1:length(Ic),uid);
            if ~isempty(cid)
                cv      = Ic(cid); 
                cv      = unique(cv);
                for ctr5=1:length(cv)
                   fc       = find(Ic==cv(ctr5));
                   cvals    = CMAT(sub2ind(size(CMAT),Ir(fc),Ic(fc)));
                   %Find largest
                   [~, ci]  = max(cvals);
                   % Delete others
                   fc(ci)   = [];
                   Ir(fc)   = [];
                   Ic(fc)   = [];
                   fl(fc)   = [];
                end
            end
            uI=unique(Ir);
            for ctr5=1:length(uI)
                idx=uI(ctr5);
                fc=find(Ir==idx);
                CLUSTS(idx,cl(idx)+1:cl(idx)+length(fc))=[Ic(fc)]';
                uC=[uC; Ic(fc)];
                cl(idx)=cl(idx)+length(fc);
                CMAT(:,Ic(fc))=0;
                nbrs=NB(CLUSTS(idx,1:cl(idx)),:);
                nbrs=nbrs(:)';
                nbrs(nbrs==0)=[];
                nbrs=setdiff(nbrs,[MaskedOut; uC]);
                inds          = CLUSTS(idx,1:cl(idx));
                FULLMAT=(Data(:,inds)'*Data(:,nbrs))./(size(Data,1)-1);
                m2=min(FULLMAT,[],1);
                CMAT(idx,:)=0;
                CMAT(idx,nbrs)=0.5+0.5*m2;
            end
        else
            RUN=0;
        end
        val=max(max(CMAT));
        ctr=ctr+length(Ic);
    end
    MAP4=zeros(size(Data,2),1);
    for ctr=1:size(CLUSTS,1);
        MAP4(CLUSTS(ctr,1:cl(ctr)))=ctr;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FAST MIN BASED   
elseif strcmp(use_dist,'fast_max')

    val=max(max(CMAT));
    RUN=1;
    ctr=0;
    CUT=0.01;
    time_run=0;
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    % relaxation parameter
    sc=0.9;
    tic
    while RUN
       % if ctr/size(Data,2)>CUT
       %     CUT=CUT+0.01;

        percent_done = ctr/size(Data,2)
        if ctr>0
            time_run     = time_run+toc/60
            approx_time_to_go   = time_run/percent_done-time_run
        end
        tic
        %end
        if val>0
            fl      = find(CMAT>=sc*val);
            [Ir Ic] = ind2sub(size(CMAT),fl);
            % See if there are several vertices that are to be assigned to different regions, if so, pick largest and delete the others 
            [u uid] = unique(Ic);
            cid     = setdiff(1:length(Ic),uid);
            if ~isempty(cid)
                cv      = Ic(cid); 
                cv      = unique(cv);
                for ctr5=1:length(cv)
                   fc       = find(Ic==cv(ctr5));
                   cvals    = CMAT(sub2ind(size(CMAT),Ir(fc),Ic(fc)));
                   %Find largest
                   [~, ci]  = max(cvals);
                   % Delete others
                   fc(ci)   = [];
                   Ir(fc)   = [];
                   Ic(fc)   = [];
                   fl(fc)   = [];
                end
            end
            uI=unique(Ir);
            for ctr5=1:length(uI)
                idx=uI(ctr5);
                fc=find(Ir==idx);
                CLUSTS(idx,cl(idx)+1:cl(idx)+length(fc))=[Ic(fc)]';
                uC=[uC; Ic(fc)];
                cl(idx)=cl(idx)+length(fc);
                CMAT(:,Ic(fc))=0;
                nbrs=NB(CLUSTS(idx,1:cl(idx)),:);
                nbrs=nbrs(:)';
                nbrs(nbrs==0)=[];
                nbrs=setdiff(nbrs,[MaskedOut; uC]);
                inds          = CLUSTS(idx,1:cl(idx));
                FULLMAT=(Data(:,inds)'*Data(:,nbrs))./(size(Data,1)-1);
                m2=max(FULLMAT,[],1);
                CMAT(idx,:)=0;
                CMAT(idx,nbrs)=0.5+0.5*m2;
            end
        else
            RUN=0;
        end
        val=max(max(CMAT));
        ctr=ctr+length(Ic);
    end
    MAP4=zeros(size(Data,2),1);
    for ctr=1:size(CLUSTS,1);
        MAP4(CLUSTS(ctr,1:cl(ctr)))=ctr;
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FAST PCA BASED   
elseif strcmp(use_dist,'fast_pca')

    val=max(max(CMAT));
    RUN=1;
    ctr=0;
    CUT=0.01;
    time_run=0;
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    % relaxation parameter
    sc=0.9;
    tic
    while RUN
       % if ctr/size(Data,2)>CUT
       %     CUT=CUT+0.01;

        percent_done = ctr/size(Data,2)
        if ctr>0
            time_run     = time_run+toc/60
            approx_time_to_go   = time_run/percent_done-time_run
        end
        tic
        %end
        if val>0
            fl      = find(CMAT>=sc*val);
            [Ir Ic] = ind2sub(size(CMAT),fl);
            % See if there are several vertices that are to be assigned to different regions, if so, pick largest and delete the others 
            [u uid] = unique(Ic);
            cid     = setdiff(1:length(Ic),uid);
            if ~isempty(cid)
                cv      = Ic(cid); 
                cv      = unique(cv);
                for ctr5=1:length(cv)
                   fc       = find(Ic==cv(ctr5));
                   cvals    = CMAT(sub2ind(size(CMAT),Ir(fc),Ic(fc)));
                   %Find largest
                   [~, ci]  = max(cvals);
                   % Delete others
                   fc(ci)   = [];
                   Ir(fc)   = [];
                   Ic(fc)   = [];
                   fl(fc)   = [];
                end
            end
            uI=unique(Ir);
            for ctr5=1:length(uI)
                idx=uI(ctr5);
                fc=find(Ir==idx);
                CLUSTS(idx,cl(idx)+1:cl(idx)+length(fc))=[Ic(fc)]';
                uC=[uC; Ic(fc)];
                cl(idx)=cl(idx)+length(fc);

                % update CMAT
                % First set column for selected index to zero
                CMAT(:,Ic(fc))=0;
                % Add elements for new neighbours of just grown region
                % fCMAT=find(CMAT(idx,:));
                nbrs=NB(CLUSTS(idx,1:cl(idx)),:);
                nbrs=nbrs(:)';
                nbrs(nbrs==0)=[];
                %nbrs=setdiff(nbrs,union(MaskedOut,union([CLUSTS(:)],fCMAT)));
                nbrs=setdiff(nbrs,[MaskedOut; uC]);

                inds          = CLUSTS(idx,1:cl(idx));
                OPTS.tol = 1e-4;
                [~, ~, D3(:,idx)] = svds(Data(:,inds)',1,'L',OPTS);
                D3(:,idx)         = D3(:,idx)-mean(D3(:,idx));
                STD               = std(D3(:,idx));
                if STD>1e-12
                    D3(:,idx)=D3(:,idx)/STD;
                else
                    error('ERROR:STD of timecourse is too small, make sure this does not happen.')
                end
                m2=(D3(:,idx)'*Data(:,nbrs))./(size(Data,1)-1);
                %distMAT(idx,nbrs)=0.5+0.5*m2;
                %CMAT(idx,nbrs)=distMAT(idx,nbrs);
                CMAT(idx,:)=0;
                CMAT(idx,nbrs)=0.5+0.5*m2;
            end
        else
            RUN=0;
        end
        val=max(max(CMAT));
        ctr=ctr+length(Ic);
    end
    MAP4=zeros(size(Data,2),1);
    for ctr=1:size(CLUSTS,1);
        MAP4(CLUSTS(ctr,1:cl(ctr)))=ctr;
    end
else
    error(['ERROR: t_grow_regions_from_seeds. Unknown method: ' use_dist]); 
end