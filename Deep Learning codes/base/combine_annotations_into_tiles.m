function [numann,percann]=combine_annotations_into_tiles(numann0,numann,percann,imlist,nblack,pthDL,outpth,sxy,dil,stile,nbg,type_emphasize,space)
if ~exist('dil','var');dil=0;end
if ~exist('nbg','var');nbg=0;end
if ~exist('stile','var');stile=10245;end
stile=stile+200;
kpall=1;


% define folder locations
outpthim=[pthDL,outpth,'im\'];
outpthlabel=[pthDL,outpth,'label\'];
outpthbg=[pthDL,outpth,'big_tiles\'];
mkdir(outpthim);mkdir(outpthlabel);mkdir(outpthbg);
imlistck=dir([outpthim,'*tif']);
nm0=length(imlistck)+1;

% create very large blank images
imH=ones([stile stile 3])*nbg;
imT=zeros([stile stile]);
nL=numel(imT);
ct=zeros([1 size(numann,2)]);
sf=sum(ct)/nL;

tic;
count=1;tcount=1;
cutoff=0.55;type0=1;
rsf=10;szzf=26;
h=zeros(szzf);h(ceil(szzf/2),ceil(szzf/2))=1;h=double(bwdist(h)<ceil(szzf/2));
while sf<cutoff
    % choose a tile containing the least prevalent class
    if rem(count,6)==1
        type=type_emphasize;
    % choose one of each class in order in a loop
    elseif rem(count+1,3)==1
        type=tcount;tcount=rem(tcount,length(ct))+1;
    else
        tmp=sum(ct,1);tmp(type0)=max(tmp);
        type=find(tmp==min(tmp),1,'first');
    end
    num=find(numann(:,type)>0);

    if isempty(num)
        numann(:,type)=numann0(:,type);
        num=find(numann(:,type)>0);
    end
    num=num(randperm(length(num),1))';

    % load annotation and mask
    imnm=imlist(num).name;
    
    pthim=[imlist(num).folder,'\'];
    pf=strfind(pthim,'\');
    pthlabel=[pthim(1:pf(end-1)),'label\'];
    
    TA=double(imread([pthlabel,imnm]));
    im=double(imread([pthim,imnm]));

    % keep only needed annotation classes
    if rem(count,3)==1;doaug=1;else;doaug=0;end % 2
    %if rem(count,3)==0;kpall=0;else;kpall=1;end
    
    [im,TA,kp]=edit_annotation_tiles(im,TA,doaug,type,ct,size(imT,1),kpall,dil);
    wsfnd=numann(num,space)>700000 & sum(numann(num,[1:space-1 space+1:end]))==0;
    numann(num,kp)=0;
    percann(num,kp)=percann(num,kp)+1; % how many times it's been used
    fx=TA~=0;
    
    if sum(fx(:))<30 || wsfnd==1
        disp('skipped');
        continue;
    end
    
    % find low density location in large tile to add annotation
    tmp=double(imT(1:rsf:end,1:rsf:end)>0);
    tmp2=imfilter(tmp,h);
    tmp=bwdist(tmp2>prctile(tmp2(:),5));
    tmp([1:10 end-9:end],:)=0;
    tmp(:,[1:10 end-9:end])=0;
    xii=find(tmp==max(tmp(:)));
    xii=xii(randperm(length(xii),1));
    [x,y]=ind2sub(size(tmp),xii);
    x=x*rsf;y=y*rsf;
    szz=size(TA)-1;
    szzA=floor(szz/2);
    szzB=szz-szzA;

    if x+szzA(1)>size(imT,2);x=x-szzA(1);end
    if y+szzA(2)>size(imT,1);y=y-szzA(2);end
    if x-szzB(1)<1;x=x+szzB(1);end
    if y-szzB(2)<1;y=y+szzB(2);end
    tmpT=imT(x-szzB(1):x+szzA(1),y-szzB(2):y+szzA(2));

    tmpT(fx)=TA(fx);
    tmpH=imH(x-szzB(1):x+szzA(1),y-szzB(2):y+szzA(2),:);
    tmpH(cat(3,fx,fx,fx))=im(cat(3,fx,fx,fx));
    imT(x-szzB(1):x+szzA(1),y-szzB(2):y+szzA(2))=tmpT;
    imH(x-szzB(1):x+szzA(1),y-szzB(2):y+szzA(2),:)=tmpH;
     
    % update total count
    if  mod(count,2)==0
        tmp=imT(101:end-100,101:end-100);
        sf=sum(tmp(:)>0)/numel(tmp);
    end 
    for p=1:size(numann,2);ct(p)=ct(p)+sum(tmpT(:)==p);end
    
    if  mod(count,150)==0 || sf>cutoff
        %figure(42);imagesc(imT);axis equal;axis off;hold on;scatter(y,x,'r*')
        tmp=histcounts(imT(:),0:size(numann,2)+1);
        ct=tmp(2:end);ct(ct==0)=1;
        %disp(round([sf*100 toc]))
    end
    
    count=count+1;type0=type;
    
end

% cut edges off tile
imH=uint8(imH(101:end-100,101:end-100,:));
imT=uint8(imT(101:end-100,101:end-100,:));
for p=1:nblack;ct(p)=sum(imT(:)==p);end
imT(imT==0)=nblack;

% save cutouts to outpth
%tmp=histcounts(imT,1:14);tmp=tmp(1:12);tmp=round(tmp./sum(tmp)*100);disp([tmp])
sz=size(imH);
for s1=1:sxy:sz(1)
    for s2=1:sxy:sz(2)
        try
            imHtmp=imH(s1:s1+sxy-1,s2:s2+sxy-1,:);
            imTtmp=imT(s1:s1+sxy-1,s2:s2+sxy-1,:);
        catch
            continue
        end
        imwrite(imHtmp,[outpthim,num2str(nm0),'.tif']);
        imwrite(imTtmp,[outpthlabel,num2str(nm0),'.tif']);
        
        nm0=nm0+1;
    end
end
nm1=dir([outpthbg,'HE*jpg']);nm1=length(nm1)+1;

% save large tiles
imwrite(imH,[outpthbg,'HE_tile_',num2str(nm1),'.jpg']);
imwrite(imT,[outpthbg,'label_tile_',num2str(nm1),'.jpg']);

end

function [im,TA,kpout]=edit_annotation_tiles(im,TA,doaug,type,ct,sT,kpall,dil)
% makes sure annotation distribution doesn't vary by more than 1%
    if doaug
        [im,TA]=augment_annotation(im,TA,1,1,1,1,0);
    else
        [im,TA]=augment_annotation(im,TA,1,1,0,0,0);
    end
    
    if kpall==0
        maxn=ct(type);
        kp=ct<=maxn*1.05;
    else
        kp=ct>=0;
    end
    % classes to always keep in the annotation tiles
    %kp(12)=1;kp(14)=1; % collagen, whitespace - whole monkey
    %kp(3)=1;kp(4)=1; % alveoli, whitespace - monkey lungs
    
    if dil
        % dilate each class a bit
        TA0=TA;
        for b=1:max(TA(:))
            tmp=TA0==b;
            tmp=imdilate(tmp,strel('disk',8));
            tmp=tmp.*(TA0==0);
            TA(tmp==1)=b;
        end
    end

    kp=[0 kp];
    tmp=kp(TA+1);

    dil=randi(15)+15;
    tmp=imdilate(tmp,strel('disk',dil));
    TA=TA.*double(tmp);
    im=im.*double(tmp);
    kpout=unique(TA);kpout=kpout(2:end);
    
    p1=min([sT size(TA,1)]);
    p2=min([sT size(TA,2)]);
    im=im(1:p1,1:p2,:);TA=TA(1:p1,1:p2);

    im=uint8(im);
    TA=uint8(TA);
end
