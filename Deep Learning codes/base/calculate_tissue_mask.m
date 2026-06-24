function [im0,TA,outpth]=calculate_tissue_mask(pth,imnm)
% creates logical image with tissue area

outpth=[pth,'TA\'];
if ~isfolder(outpth);mkdir(outpth);end

try im0=imread([pth,imnm,'.tif']);
catch
    im0=imread([pth,imnm,'.jpg']);
end
if exist([outpth,imnm,'.tif'],'file')
    TA=imread([outpth,imnm,'.tif']);
    return;
end

disp('      calculating TA image')
if exist([outpth,'TA_cutoff.mat'],'file')
    load([outpth,'TA_cutoff.mat'],'cts')
    ct=mean(cts);
else
    ct=210;
end

TA=im0(:,:,2)<ct;
TA=imclose(TA,strel('disk',1));
TA=bwareaopen(TA,10);
TA=~bwareaopen(~TA,10);
imwrite(TA,[outpth,imnm,'.tif']);

