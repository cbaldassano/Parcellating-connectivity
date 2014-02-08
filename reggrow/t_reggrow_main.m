function [CLUST, TREE, REG_INDS, BORDERS, ROI_var, SMOOTH, OUT, REGS] = t_reggrow_main(DATA, NBRS, VERTS, NC);
% [CLUST, TREE, REG_INDS, BORDERS, ROI_var, SMOOTH, OUT, REGS] = t_reggrow_main(DATA, NBRS, VERTS, NC);
%
% spatially constrained brain parcellation
% DATA - dx by dt data matrix, where each row is a time series associated
%        with a location
% NBRS - dx by dn neighbourhood structure. Each location in DATA (i.e. each row in DATA) has a row in NBRS
%        The numbers in each row in NBRS are the indices of the rows in DATA that are spatail neighbours 
% VERTS - dx by 3 matrix of locations. Specifies the spatial coordinates in
%         xyz for each time series in DATA
% NC    - Number of clusters to return
%

%
filtwidth=1;

%% FIND NEIGHBORHOOD STRUCTURE

[dx dt]=size(DATA);

%% FIND ROI STRUCTURE

%display(['Defining ROIs of radius 1mm, 2mm, 3mm'])
ROI1mm  = zeros(size(NBRS));
ROI2mm  = zeros(size(NBRS));
ROI3mm  = zeros(size(NBRS));
%ROI9mm  = zeros(size(NBRS));
nr1   = size(ROI1mm,2);
nr2   = size(ROI2mm,2);
nr3   = size(ROI3mm,2);
%nr9   = size(ROI9mm,2);
nroi1 = zeros(size(ROI1mm,1),1);
nroi2 = zeros(size(ROI2mm,1),1);
nroi3 = zeros(size(ROI3mm,1),1);
%nroi9 = zeros(size(ROI9mm,1),1);
tic
for ctr=1:size(ROI3mm,1)
    if rem(ctr,1000)==0
        %display(['ROI definition iteration: ' num2str(ctr)])
    end
    closeverts1 = [];
    closeverts2 = [];
    closeverts3 = [];
   % closeverts9 = [];
    nbrs = NBRS(ctr,:);
    nbrs(nbrs==0)=[];
    nbrs    = union(nbrs, ctr);
    CHECKED = nbrs;
    centre  = VERTS(ctr,:);        
    nbrdist = sqrt(sum([ (VERTS(nbrs,1)-centre(1)).^2 (VERTS(nbrs,2)-centre(2)).^2 (VERTS(nbrs,3)-centre(3)).^2],2));
    newnbr1  = nbrs(nbrdist<=1);
    newnbr2  = nbrs(nbrdist<=2);
    newnbr3  = nbrs(nbrdist<=3);
   % newnbr9  = nbrs(nbrdist<=9);
    while ~isempty(newnbr3)
        closeverts1 = [closeverts1 newnbr1];
        closeverts2 = [closeverts2 newnbr2];
        closeverts3 = [closeverts3 newnbr3];
    %    closeverts9 = [closeverts9 newnbr9];
        nbrs = NBRS(nbrs,:);
        nbrs=nbrs(:);
        nbrs(nbrs==0) = [];
        nbrs = setdiff(nbrs(:)',CHECKED);
        CHECKED = [CHECKED nbrs];
        nbrdist = sqrt(sum([ (VERTS(nbrs,1)-centre(1)).^2 (VERTS(nbrs,2)-centre(2)).^2 (VERTS(nbrs,3)-centre(3)).^2],2));
        newnbr1  = nbrs(nbrdist<=1);
        newnbr2  = nbrs(nbrdist<=2);
        newnbr3  = nbrs(nbrdist<=3);
    %    newnbr9  = nbrs(nbrdist<=9);
    end
    nroi1 (ctr) = length(closeverts1);
    nroi2 (ctr) = length(closeverts2);
    nroi3 (ctr) = length(closeverts3);
   % nroi9 (ctr) = length(closeverts9);
    if nr1<nroi1(ctr) ;
        ROI1mm=[ROI1mm zeros(size(ROI1mm,1),nroi1(ctr)-nr1)];
        nr1=nroi1(ctr) ;
    end
        ROI1mm(ctr,1:nroi1(ctr))=closeverts1;
    if nr2<nroi2(ctr) ;
        ROI2mm=[ROI2mm zeros(size(ROI2mm,1),nroi2(ctr)-nr2)];
        nr2=nroi2(ctr) ;
    end
        ROI2mm(ctr,1:nroi2(ctr))=closeverts2;
    if nr3<nroi3(ctr) ;
        ROI3mm=[ROI3mm zeros(size(ROI3mm,1),nroi3(ctr)-nr3)];
        nr3=nroi3(ctr) ;
    end
        ROI3mm(ctr,1:nroi3(ctr))=closeverts3;
    %if nr9<nroi9(ctr) ;
    %    ROI9mm=[ROI9mm zeros(size(ROI9mm,1),nroi9(ctr)-nr9)];
    %    nr9=nroi9(ctr) ;
    %end
    %    ROI9mm(ctr,1:nroi9(ctr))=closeverts9;
end
%toc
ROI1mm = [ROI1mm nroi1]; % roi = ROI(c,1:ROI(c,end)); returns roi of vertex c
ROI2mm = [ROI2mm nroi2];
ROI3mm = [ROI3mm nroi3];
%ROI9mm = [ROI9mm nroi9];

%% Define convolution kernal


%display('Creating convolution kernal matrix')
% Create a matrix with filter weights. Once this is done, filter can be
% applied repeatedly.
FILT_WEIGHTS=sparse([],[],[],dx,dx, prod(size(ROI3mm)));
for ctr=1:dx
    if rem(ctr,1000)==0
        display(['Iteration: ' num2str(ctr)])
    end
    c        = ctr;
    roi      = ROI3mm(c,1:ROI3mm(c,end));
  % roi      = intersect(roi,fM);
    centre   = VERTS(c,:);        
    nbrdist  = sqrt(sum([ (VERTS(roi,1)-centre(1)).^2 (VERTS(roi,2)-centre(2)).^2 (VERTS(roi,3)-centre(3)).^2],2));
    W        = exp(-nbrdist.^2/(2*filtwidth^2));
    W        = W./sum(W);
    FILT_WEIGHTS(roi,c)=W;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DOWNSAMPLE

[DATA_DS]  = t_downsample(squeeze(DATA), 0);
VoxSize_DS = size(DATA_DS);
clear data

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PREPROCESSING
DEMEAN=1;
VoxVarNorm=1;
[Data_PRE, Mask] = t_prepro(DATA_DS, DEMEAN, VoxVarNorm);
fM=find(Mask);
Data=zeros(dt,dx);
Data(:,fM)=Data_PRE;

clear Data_PRE
%display('done')

fM=find(Mask);
Data=Data/norm(Data(:,fM(1)));

%% Calculate ROI variance

%display(['Estimating ROIs variation in radius of 3mm'])
ROI_var=zeros(size(ROI3mm,1),1);
for ctr=1:length(fM)
    if rem(ctr,1000)==0
        %display(['ROI variation iteration: ' num2str(ctr)])
        %figure(1);plot(ROI_var(1:c))
    end
    c   = fM(ctr);
    roi = ROI3mm(c,1:ROI3mm(c,end));
    roi = intersect(roi,fM);
    %roi(roi==c)=[];
    D   = Data(:,roi);
    mD  = mean(D,2);
    ERR=D-repmat(mD,[1 size(D,2)]);
    ROI_var(c) = var(ERR(:));
end

%% SMOOTH with Gaussian kernal
% Spatial smoothing with gaussina kernal
% Nadaraya-Watson type kernal. 
% Method is similar to: Moo K. Chung similar to Moo K. Chung, but using euclidean space distances which,
% for close points, should be a good approximation ot geodesic distances. Better 
% approximation might be distance along edges, whilst even better might be an aver
% age between edge distance and euclidean distance.

%%
%display(['Averaging (spatial) data with Gaussian kernal of width 3mm'])
IN=ROI_var';
fNM=find(Mask==0);
IN(fNM)=max(IN);

SM=1;
for ctr=1:SM
    IN = IN*FILT_WEIGHTS;
    IN(fNM)=max(IN);
end
SMOOTH=IN;

%% FIND LOCAL MINIMA

[OUT, OUT2]=t_find_local_min(SMOOTH', NBRS, 0);
OUT(fNM)=0;

while length(unique(OUT))>5000
    IN = IN*FILT_WEIGHTS;
    IN(fNM)=max(IN);
    SMOOTH=IN;
    [OUT, OUT2]=t_find_local_min(SMOOTH', NBRS, 0);
    OUT(fNM)=0;
end


%%
%MAP=zeros(size(Data,2),1);

uO=unique(OUT);
uO(uO==0)=[];
D=zeros(size(Data,1),length(uO));
ctr=1;
for ctr1=uO'
   c=find(OUT==ctr1);
   roi = ROI3mm(c,1:ROI3mm(c,end));
   %MAP(roi)=ctr;
   D(:,ctr)=mean(Data(:,roi),2);
   D(:,ctr)=D(:,ctr)-mean(D(:,ctr));
   D(:,ctr)=D(:,ctr)/norm(D(:,ctr));
   OUT(c)=ctr;
   ctr=ctr+1;
end



%% STEP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Average bordermaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%    


OUT  = OUT  .* Mask;
fOUT=find(OUT);
old_OUT2_CENTRES=OUT2(fOUT);
OUT2 = OUT2 .* Mask;
new_OUT2_CENTRES=OUT2(fOUT);
fCHANGE=find(old_OUT2_CENTRES-new_OUT2_CENTRES);
nums=unique(old_OUT2_CENTRES(fCHANGE));
for ctr=1:length(nums)
 if isempty(find(new_OUT2_CENTRES==nums(ctr), 1))
     OUT2(OUT2==nums(ctr))=0;
 end
end
OUT(fOUT)=OUT2(fOUT); % Make sure we use the same labels in both, just in case.
SL=setdiff(unique(OUT2) ,unique(OUT));
for ctr=1:length(SL)
    OUT2(OUT2==SL(ctr))=0; 
end

uO=unique(OUT);
for ctr=uO(1):length(uO)-1
    OUT(OUT==uO(ctr+1))=ctr;
    OUT2(OUT2==uO(ctr+1))=ctr;
end

if length(unique(OUT))~=length(unique(OUT2))
    error('ERROR: OUT and OUT2 should have the same cluster labels! Something is wrong with the local minima assignment and labelling!')
end
if sum(unique(OUT)-unique(OUT2))
    error('ERROR: OUT and OUT2 should have the same cluster labels! Something is wrong with the local minima assignment and labelling!')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% GROW SEED REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%


REGS=t_grow_regions_from_seeds(Data, OUT, OUT2, Mask, NBRS, 1 ,'fast');

MAPt=REGS;

uO=unique(OUT);
if uO(1)==0;
    uO(1)=[];
end
for ctr=uO';
    fO=find(OUT==ctr);
    MM=REGS(fO);
    for ctr1=1:length(MM)
        fM=find(REGS==MM(ctr1));
        MAPt(fM)=ctr;
    end
end
REGS=MAPt;
%%      

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% EXTRACT SEED TIMECOURSES %%%%%%%%%%%%%%%%%%%%%%%%%
    
inds=unique(OUT2);
if inds(1)==0
    inds=inds(2:end)';
end
D=zeros(size(Data,1),length(inds));
ctr=1;
NOTIN=[];
for ctr1=inds    
    fi=find(OUT2(:,1)==ctr1);
    if ~isempty(fi)
        D(:,ctr)=mean(Data(:,fi),2);
        D(:,ctr)=D(:,ctr)-mean(D(:,ctr));
        STD=std(D(:,ctr));
        if STD<=1e-12 %Exclusde time-course with very small std THIS SHOULD BE NORMALISED IN SOME WAY TO GUARANTEE ROBUSTNES!
            NOTIN=[NOTIN ctr];
        else
            D(:,ctr)=D(:,ctr)/STD;
        end
        ctr=ctr+1;
    end
end
D(:,NOTIN)=[];
REG_INDS=inds;
REG_INDS(NOTIN)=[];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PRE-CLUSTER SEEDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Generate local neighbourhood structure %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

maxnb=5;
regNB=zeros(length(REG_INDS),5);
uM=unique(REGS);
if uM(1)==0;
    uM(1)=[];
end
for ctr=1:length(uM)
    ff=find(REGS(:,1)==uM(ctr));
    nbrs=NBRS(ff,:);
    nbrs=nbrs(:);
    nbrs=unique(nbrs);
    nbrs(nbrs==0)=[];
    nbrs=setdiff(nbrs,ff);
    tmpnb=unique(REGS(nbrs));
    tmpnb(tmpnb==0)=[];
    ltn=length(tmpnb);
    if ltn>maxnb
        regNB=[regNB zeros(length(REG_INDS), (ltn-maxnb))];
        maxnb=ltn;
    end
    regNB(ctr,1:ltn)=tmpnb;
end

%% Cluster

DD=D'*D/(size(Data,1)-1);
dist=(2-2*DD).^0.5; % Use this so that dist is interpretable as euclidean
%dist=(2-2*DD);
dist=dist-diag(diag(dist));
tclust='local_ward';
TREE = t_local_linkage(squareform(dist),regNB,'wa');

lNC=length(NC);
CLUST=zeros(dx,lNC);
BORDERS=zeros(dx,lNC);
for nc=1:lNC
    num_clust=NC(nc);
    C = cluster(TREE,'maxclust',num_clust);

    %% Make cluster map

    MAP6=zeros(size(Data,2),1);
    for ctr=1:length(REG_INDS)
        ff=find(REGS(:,1)==REG_INDS(ctr));
        MAP6(ff,1)=C(ctr);
    end
   % ciftiout1.cdata=zeros(size(ciftiout1.cdata,1),1)
    % Randomly permute labels for better visualisation
    IN=MAP6;
    uO=unique(IN);
    uO(uO==0)=[];

    RP=randperm(length(uO));
    for ctr=1:length(uO)
        f=find(IN==uO(ctr));
        CLUST(f,nc)=RP(ctr);
    end

    for ctr=1:size(BORDERS,1);
        nbrs=NBRS(ctr,:);
        nbrs(nbrs==0)=[];
        c=MAP6(ctr);
        cc=MAP6(nbrs);
        cd=setdiff(cc,c);
        if ~isempty(cd)
           BORDERS(ctr,nc)=1; 
        end
    end
end