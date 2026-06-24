function outpth=classify_images(pthim,pthDL,colHE,col)
% pthim: path to tif images to classify
% pthDL: location where net.mat file with model is saved
% colHE (optional): save classification superimposed on H&E. default = yes
% col (optional): save colorized classification. default = no
if ~exist('colHE','var');colHE=1;end
if ~exist('col','var');col=0;end

load([pthDL,'net.mat'],'net','nm','cmap','sxy','nblack','nwhite');
outpth=[pthim,'classification_',nm,'/'];
mkdir(outpth);

b=100;
imlist=dir([pthim,'*tif']);
if isempty(imlist);imlist=dir([pthim,'*jp2']);end
if isempty(imlist);imlist=dir([pthim,'*jpg']);end
disp(' ')

x=tic;
for kk=1:length(imlist)
    tic;nm=imlist(kk).name;
    disp(['  Starting classfication of image ',num2str(kk),' of ',num2str(length(imlist)),': ',nm])
    if exist([outpth,nm(1:end-3),'tif'],'file')
        disp('    image already classified with this model');
        continue;
    end
    try
    im=imread([pthim,nm]);
    catch 
        continue
    end

    try
        TA=imread([pthim,'TA/',nm,'tif']);
        imfill(TA,'holes');
    catch
        TA=rgb2gray(im)<220;
        imfill(TA,'holes');
    end
    
    % pad image so we classify all the way to the edge
    im=padarray(im,[sxy+b sxy+b],0,'both');
    TA=padarray(TA,[sxy+b sxy+b],1,'both');
    
    imclassify=zeros(size(TA));
    sz=size(im);
    for s1=1:sxy-b*2:sz(1)-sxy
        for s2=1:sxy-b*2:sz(2)-sxy
            tileHE=im(s1:s1+sxy-1,s2:s2+sxy-1,:);
            tileTA=TA(s1:s1+sxy-1,s2:s2+sxy-1,:);
            if sum(tileTA(:))<100
                tileclassify=zeros(size(tileTA));
            else
                tileclassify=semanticseg(tileHE,net);
            end
            tileclassify=tileclassify(b+1:end-b,b+1:end-b,:);
            imclassify(s1+b:s1+sxy-b-1,s2+b:s2+sxy-b-1)=tileclassify;

%             tileclassify=double(tileclassify);
%             tileclassify(tileclassify==nblack | tileclassify==0)=nwhite;
%             imcolor=uint8(cat(3,am(tileclassify),bm(tileclassify),cm(tileclassify)));
%             tileHE=tileHE(b+1:end-b,b+1:end-b,:);
%             figure(18),
%                 subplot(1,2,1),imshow(tileHE)
%                 subplot(1,2,2),imshow(imcolor)
        end
    end
   
    % remove padding
    im=im(sxy+b+1:end-sxy-b,sxy+b+1:end-sxy-b,:);
    imclassify=imclassify(sxy+b+1:end-sxy-b,sxy+b+1:end-sxy-b,:);
    
    disp([kk length(imlist) round(toc)])
    imclassify(imclassify==nblack | imclassify==0)=nwhite; % make black class and zeros class whitespace
    imwrite(uint8(imclassify),[outpth,nm(1:end-3),'tif']);
    
    % make colorized image
    if col==1
        outpthcolor=[outpth,'color/'];
        if ~isfolder(outpthcolor);mkdir(outpthcolor);end
        am=cmap(:,1);bm=cmap(:,2);cm=cmap(:,3);
        imcolor=uint8(cat(3,am(imclassify),bm(imclassify),cm(imclassify)));
        imwrite(imcolor,[outpthcolor,nm(1:end-3),'tif']);
    end

    % make color image overlayed on H&E
    if colHE==1
        outpth2=[outpth,'check_classification/'];
        if ~isfolder(outpth2);mkdir(outpth2);end
        make_check_annotation_image(im,imclassify,cmap,2,[outpth2,nm(1:end-3),'jpg']);
    end

    % display the first image in the series
    if kk==1 && ~isempty(cmap)
        am=cmap(:,1);bm=cmap(:,2);cm=cmap(:,3);
        imcolor=uint8(cat(3,am(imclassify),bm(imclassify),cm(imclassify)));
        figure;subplot(1,2,1);imshow(im);subplot(1,2,2);imshow(imcolor)
        ha=get(gcf,'children');linkaxes(ha);
        pause(0.2)
    end

end
disp(['  Total time for classification: ',num2str(round(toc(x)/60)),' minutes'])

