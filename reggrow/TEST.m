% Random reggrow example

dx=100;
dt=10;

addpath([pwd '/hlpfunc'])

%% Generate random data

DATA=randn(dx,dt);

%% Generate some neighbourhood structure, this must be a matrix, zero fill rows if there are fewer neighbours
N=ceil(sqrt(dx));
NBRS=zeros(dx,4);

for ctr=1:dx
   [col row]=ind2sub([N N],ctr);
   try
       i1=sub2ind([N N],col-1,row);
   catch
       i1=0;
   end
   try
       i2=sub2ind([N N],col,row-1);
   catch
       i2=0;
   end
   try
       i3=sub2ind([N N],col+1,row);
   catch
       i3=0;
   end
   try
       i4=sub2ind([N N],col,row+1);
   catch
       i4=0;
   end
   nbrs=[i1 i2 i3 i4];
   nbrs(nbrs<0)=0;
   fn=find(nbrs>0);
   NBRS(ctr,1:length(fn))=sort(nbrs(fn));
end
NBRS(NBRS>dx)=0;
%% Generate random vertex locations
[cols rows]=ind2sub([N N],1:dx)
VERTS=[cols' rows' rand(dx,1)];

%% RUN REGGROW

NC=10; % Number of clusters
[CLUST, TREE, REG_INDS, BORDERS, ROI_var, SMOOTH, OUT, REGS] = t_reggrow_main(DATA, NBRS, VERTS, NC);

%% 

figure(1)
imagesc(reshape(CLUST,N,N))