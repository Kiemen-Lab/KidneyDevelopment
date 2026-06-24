function validate_annotations_highres(pth,nm,I0,J0,imnm,classNames,cmap2,classcheck)
% pth = path to xml files
% pthim = path to H&E tif images
% classcheck = CNN class you want to make bounding boxes for
% classNames = title of each CNN class from deeplab code
% cmap = color of each class

if classcheck==0
    disp('    skipping validation images')
    return;
else
    disp(['    3. of 4. Saving validation images for ',num2str(length(classcheck)),' classes'])
end

outpth=[pth,'validate_each_class_',nm,'/'];
% if exist([outpth,imnm,'_validate.mat'],'file');return;end

% ds=2;d=1; 
ds=1;d=1;
I0=im2double(I0(1:ds:end,1:ds:end,:));
for tp=classcheck
    J=double(J0(1:ds:end,1:ds:end));
    
    % set up names and folders
    nmcheck=char(classNames(tp));
    outim=[outpth,'check_',nmcheck,'/'];
    if ~isfolder(outim);mkdir(outim);end

    % save bounding boxes containing tissue of interest
    % tmp=imclose(J==tp,strel('disk',5));
    tmp =  J==tp;
    tmp=imfill(tmp,'holes');
    tmp=bwareaopen(tmp,100);
    L=bwlabel(tmp);
    b=unique(L(:));
    disp(['      saving ',num2str(length(b)-1),' ',nmcheck,' annotations for validation'])
    if length(b)==1;continue;end

    % make color mask
    cm = zeros(size(cmap2));
    cm(2:end,1) = 0.75;
    cm(tp+1,:)=[1 1 1];
    J1=cm(J+1,1);J1=reshape(J1,size(J));
    J2=cm(J+1,2);J2=reshape(J2,size(J));
    J3=cm(J+1,3);J3=reshape(J3,size(J));
    mask=cat(3,J1,J2,J3);
    I=im2uint8((double(I0)*0.75)+(double(mask)*0.25));
    clearvars J1 J2 J3

    pause(2)
    for pk=b(2:end)'
        tmp=double(L==pk);
        a=sum(tmp,1);b=sum(tmp,2);
        rect=[find(a,1,'first') find(a,1,'last') find(b,1,'first') find(b,1,'last')];
        rect=[max([rect(1)-20 1]) min([rect(2)+20 size(tmp,2)])...
            max([rect(3)-20 1])  min([rect(4)+20 size(tmp,1)])];
        % make label and  H&E bounding boxes
        tmp=J(rect(3):rect(4),rect(1):rect(2));
        if sum(tmp(:))==0;break;end
        tmpim=I(rect(3):rect(4),rect(1):rect(2),:);
        imwrite(uint8(tmpim),[outim,imnm,'_',num2str(pk,'%04.f'),'.jpg']);
        pause(0.05)
    end

end

% save a random variable when all boxes are saved for this image so we can skip the image next time
save([outpth,imnm,'_validate.mat'],'d');

