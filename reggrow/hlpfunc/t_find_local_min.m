function [OUT, OUT2]=t_find_local_min(dat, NB, JOIN)

[a b]=size(dat);
TRANS=0;
if a<b
    dat=dat';
    [a b]=size(dat);
    TRANS=1;
end

OUT =zeros(size(dat));
OUT2=zeros(size(dat));


for ctr1=1:b
    ctr3=1;
    ctr4=1;
    for ctr2=1:a
        nbrs=NB(ctr2,:);
        nbrs(nbrs==0)=[];
        min_n=min(dat(nbrs,ctr1));
        if dat(ctr2,ctr1) <= min_n
            OUT(ctr2,ctr1)=ctr3;
            mnbrs=max(OUT2(nbrs,ctr1)); %CB changed OUT->OUT2
            if mnbrs==0
                OUT2(ctr2,ctr1)=ctr4;
                OUT2(nbrs,ctr1)=ctr4;
                ctr4=ctr4+1;
            else
                OUT2(ctr2,ctr1)=mnbrs;
                OUT2(nbrs,ctr1)=mnbrs;
            end
            ctr3=ctr3+1;
        end
    end
    % Join touching OUT2 clusters
    if JOIN
        for ctr5=1:max(OUT2(:,ctr1))
           fc=find(OUT2(:,ctr1)==ctr5);
           nbrs=NB(fc,:);
           nbrs=nbrs(:);
           %nbrs(nbrs==0)=[];
           nbrs=setdiff(nbrs,union(fc,0));
           fn=find(OUT2(nbrs,ctr1)~=0);
           for ctr6=1:length(fn)
               nval=OUT2(nbrs(fn(ctr6)),ctr1);
               if nval==0
                   display('ERROR!')
               else
                   OUT2(OUT2(:,ctr1)==nval)=ctr5;
               end
           end
        end
    end
%     for ctr5=1:max(OUT2(:,ctr1))
%        fc=find(OUT2(:,ctr1)==ctr5);
%        nbrs=NB(fc,:);
%        nbrs=nbrs(:);
%        nbrs(nbrs==0)=[];
%        fn=find(OUT2(nbrs,ctr1)~=0);
%        for ctr6=1:fn
%            nval=OUT2(nbrs(fn(ctr6)),ctr1);
%            if nval==0
%                display('ERROR!')
%            else
%                OUT2(OUT2(:,ctr1)==nval)=ctr5;
%            end
%        end
%     end
end

if TRANS
    OUT=OUT';
    OUT2=OUT2';
end