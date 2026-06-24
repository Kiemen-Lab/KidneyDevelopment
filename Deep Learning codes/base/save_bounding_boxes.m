function [numann,ctlist]=save_bounding_boxes(im,outpth,nm0,numclass)
disp('    4. of 4. Creating bounding box tiles of all annotations')

im=double(im);
try
    imlabel=double(imread([outpth,'view_annotations.tif']));
catch
    imlabel=double(imread([outpth,'view_annotations_raw.tif']));
end


pthbb=[outpth,nm0,'_boundbox/'];
pthim=[pthbb,'im/'];
pthlabel=[pthbb,'label/'];
if isfolder('pthim');rmdir(pthim);rmdir(pthlabel);end
mkdir(pthim);mkdir(pthlabel);

tmp=imclose(imlabel>0,strel('disk',10));
tmp=imfill(tmp,'holes');
tmp=bwareaopen(tmp,300);
L=bwlabel(tmp);
numann=zeros([max(L(:)) numclass]);
for pk=1:max(L(:))
    tmp=double(L==pk);
    a=sum(tmp,1);b=sum(tmp,2);
    rect=[find(a,1,'first') find(a,1,'last') find(b,1,'first') find(b,1,'last')];
    tmp=tmp(rect(3):rect(4),rect(1):rect(2));
    
    % make label and  H&E bounding boxes
    tmplabel=imlabel(rect(3):rect(4),rect(1):rect(2)).*tmp;
    tmpim=im(rect(3):rect(4),rect(1):rect(2),:);
    nm=num2str(pk,'%05.f');
    imwrite(uint8(tmpim),[pthim,nm,'.tif']);
    imwrite(uint8(tmplabel),[pthlabel,nm,'.tif']);

    for anns=1:numclass
       numann(pk,anns)=sum(tmplabel(:)==anns);           
    end
end
ctlist=dir([pthim,'*tif']);

bb=1; % indicate that xml file is fully analyzed
if exist([outpth,'annotations.mat'],'file')
    save([outpth,'annotations.mat'],'numann','ctlist','bb','-append');
else
    save([outpth,'annotations.mat'],'numann','ctlist','bb');
end
