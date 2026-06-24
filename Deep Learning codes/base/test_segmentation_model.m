function test_segmentation_model(pthDL,pthtest,pthtestim)
disp('Testing the accuracy of the model...')
load([pthDL,'net.mat'],'nblack','nwhite','classNames');

pthtestdata=[pthtest,'data/'];
load_annotation_data(pthDL,pthtest,pthtestim);
pthclassifytest=classify_images(pthtestim,pthDL);

classNames=classNames(1:end-1);
numclass=nblack-1;

num_xml = dir([pthtest,'*xml']);
plist=dir(pthtestdata);
 if length(num_xml)~= length(plist)-2
    disp('Verify data folder');

else
    pDL=[];ptrue=[];
    for k=1:length(plist)
        tic;
        pth=[pthtestdata,plist(k).name,'/'];
        if exist([pth,'view_annotations.tif'],'file') || exist([pth,'view_annotations_raw.tif'],'file')
            try
                J0=double(imread([pth,'view_annotations.tif']));
            catch
                J0=double(imread([pth,'view_annotations_raw.tif']));
            end
            imDL=imread([pthclassifytest,plist(k).name,'.tif']);
            
            % remove small pixels
            for b=1:max(J0(:))
                tmp=J0==b;
                J0(J0==b)=0;
                tmp=bwareaopen(tmp,25);
                J0(tmp==1)=b;
            end
            
            % get true and predicted class at testing annotation locations
            L=find(J0>0);
            ptrue=cat(1,ptrue,J0(L));
            pDL=cat(1,pDL,imDL(L));
        end
    end
    pDL(pDL==nblack)=nwhite;
    
    % normalize to the minimum number of pixels, rounded to nearest 1000
    count_ann=histcounts(ptrue,numclass);
    comp_ann=(count_ann/max(count_ann)*100);
    km=min(count_ann);
    
    % display number of pixels of each class in testing
    disp('Calculating total number of pixels in the testing dataset...')
    for b=1:length(count_ann)
        if comp_ann(b)==100
            disp(['  There are ',num2str(count_ann(b)),' pixels of ',char(classNames(b)),...
                '. This is the most common class.'])
        else
            disp(['  There are ',num2str(count_ann(b)),' pixels of ',char(classNames(b)),...
                ', ',num2str(ceil(comp_ann(b))),'% of the most common class.'])
        end
    end
    
    if ~isempty(find(count_ann==0,1))
        for bb=find(count_ann==0)
            disp(' ');disp([' No testing annotations exist for class ',char(classNames(bb)),'.'])
        end
        error(' Cannot make confusion matrix. Please add testing annotations of missing class(es).')
    end
    
    if km<100
        km=floor(km/10)*10;
    elseif km<1000
        km=floor(km/100)*100;
    else
        km=floor(km/1000)*1000;
    end
    
    for bb=find(count_ann<15000)
        disp(' ')
        disp(['  Only ',num2str(count_ann(bb)),' testing pixels of ',char(classNames(bb)),' found.'])
        disp('    We suggest a minimum of 15,000 pixels for a good assessment of model accuracy.')
        disp('    Confusion matrix may be misleading.')
    end
    
    ptrue2=[];
    pDL2=[];
    for k=unique(ptrue)'
        % if k==9 | k==11
        %     continue
        % end
        a=find(ptrue==k);
        b=randperm(length(a),km);
        ptrue2=[ptrue2;ptrue(a(b))];
        pDL2=[pDL2;pDL(a(b))];
    end
    
    % confusion matrix with equal number of pixels of each class
    Dn=zeros([max(ptrue) max(pDL)]);
    for a=1:max(ptrue2)
        for b=1:max(pDL2)
        %      if a==9 | b==11 | b==9 |a==11
        %     continue
        % end
            tmp1=ptrue2==a;
            tmp2=pDL2==b;
            Dn(a,b)=sum(tmp1 & tmp2);
        end
    end
    Dn(isnan(Dn))=0;
    plot_confusion_matrix(Dn,classNames,pthDL)
end



