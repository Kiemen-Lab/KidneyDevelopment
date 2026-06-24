function xyc=combine_z_projections(vol,nums,contrasts,cmap0,tt,titles)
if ~exist('cmap0','var')
    cmap0=[121 248 252;...  % 1 islet
        40    40    255;... % 2 duct
        80 237 80;...       % 3 blood vessel
        255  255  0;...     % 4 fat
        149 35  184;...     % 5 acinus
        255 194 245;...     % 6 connective tissue
        255 255 255;...     % 7 whitespace
        255  0  0;...       % 8 PanIN
        240 159 10;...      % 9 PDAC 
        80 65 35;...        % 10 nerve
        125 125 125];   % 11 lymph
end
cmap0=cmap0/255;
% 
xyc=zeros([size(vol(:,:,1)) 3]);


for k=nums
    % make colormap
    cmap=cmap0(k,:);
    a=linspace(0,cmap(1),100)';
    b=linspace(0,cmap(2),100)';
    c=linspace(0,cmap(3),100)';
    C=[a b c]*contrasts(k); % enhance contrast
    C(C>1)=1;
    
    % make rgb z-projection
    xy0=squeeze(sum(vol==k,3)); 
    xy0=round(xy0./max(xy0(:))*99)+1;
    xy0(isnan(xy0))=1;
    xy=cat(3,C(xy0,1),C(xy0,2),C(xy0,3));
    xy=reshape(xy,[size(xy0) 3]);
    xyc=xyc+xy;
end




figure;imshow(xyc);
if exist('tt','var');title(tt);end
if exist('titles','var')
    cmap2=cmap0(nums,:)*255;
    titles2=titles(nums);
    make_cmap_legend(cmap2,titles2)
end