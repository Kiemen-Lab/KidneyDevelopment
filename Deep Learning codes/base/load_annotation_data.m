function [ctlist0,numann0]=load_annotation_data(pthDL,pth,pthim,classcheck)
    if ~exist('classcheck','var');classcheck=0;end
    if isempty(classcheck);classcheck=0;end
    disp(' ');
    disp('Importing annotation data...')
    load([pthDL,'net.mat'],'WS','umpix','cmap','classNames','nm','nwhite');
    cmap2=cat(1,[0 0 0],cmap)/255;
    numclass=max(WS{3});
    
    imlist=dir([pth,'*.xml']);
    numann0=[];ctlist0=[];
    outim=[pth,'check_annotations\'];mkdir(outim);
    
    % first, check that all images exist
    % for kk=1:length(imlist)
    %     imnm=imlist(kk).name(1:end-4);
    %     exx=exist([pthim,imnm,'.tif'],'file') || exist([pthim,imnm,'.jpg'],'file');
    %     if exx==0
    %         error([' Cannot find a tif or jpg file for xml file: ',imnm,'.xml'])
    %     end
    % end

    % for each annotation file
    for kk=1:length(imlist)
        tic;
        % set up names
        imnm=imlist(kk).name(1:end-4);
        
        disp(['  Image ',num2str(kk),' of ',num2str(length(imlist)),': ',imnm])
        outpth=[pth,'data\',imnm,'\'];
        matfile=[outpth,'annotations.mat'];

        % check if model parameters have changed
        reload_xml=check_if_model_parameters_changed(matfile,WS,umpix,nwhite,pthim);
        % reload_xml =1;
        % skip if file hasn't been updated since last load
        dm='';bb=0;date_modified=imlist(kk).date;
        if exist(matfile,'file');load(matfile,'dm','bb');end
        if contains(dm,date_modified) && bb==1 && reload_xml==0
            disp('    annotation data previously loaded')
            load([outpth,'annotations.mat'],'numann','ctlist');
            numann0=cat(1,numann0,numann);
            ctlist0=cat(1,ctlist0,ctlist);
            continue;
        end
        % 
        % if ~contains(imlist(kk).name,{'34954_278','34957','34958'}) 
        %     disp('    annotation data previously loaded')
        %     load([outpth,'annotations.mat'],'numann','ctlist');
        %     numann0=cat(1,numann0,numann);
        %     ctlist0=cat(1,ctlist0,ctlist);
        %     continue;
        % end
     % if ~contains(imlist(kk).name,{'34954_124_126'}) 
     %            disp('    annotation data previously loaded')
     %            load([outpth,'annotations.mat'],'numann','ctlist');
     %            numann0=cat(1,numann0,numann);
     %            ctlist0=cat(1,ctlist0,ctlist);
     %            continue;
     %        end

        if isfolder(outpth);delete([outpth,'/*']);end
        mkdir(outpth);
        
        % 1 read xml annotation files and saves as mat files
        import_xml(matfile,[pth,imnm,'.xml'],date_modified);
        save(matfile,'WS','umpix','nwhite','pthim','-append');
        
         % 2 fill annotation outlines and delete unwanted pixels
        [I0,TA,~]=calculate_tissue_mask(pthim,imnm);
        J0=save_annotation_mask(I0,outpth,WS,umpix,TA,1); 
        imwrite(uint8(J0),[outpth,'view_annotations.tif']);
        
        % save validation images
        validate_annotations_highres(pth,nm,I0,J0,imnm,classNames,cmap2,classcheck)
        
        % show mask in color: J2=J+1;cmap3=cmap2(min(J2(:)):max(J2(:)),:);figure,imshow(J2,cmap3);
        I=im2double(I0(1:2:end,1:2:end,:));J=double(J0(1:2:end,1:2:end,:));
        J1=cmap2(J+1,1);J1=reshape(J1,size(J));
        J2=cmap2(J+1,2);J2=reshape(J2,size(J));
        J3=cmap2(J+1,3);J3=reshape(J3,size(J));
        mask=cat(3,J1,J2,J3);
        I=(I*0.5)+(mask*0.5);
        imwrite(im2uint8(I),[outim,imnm,'.jpg']);
        clearvars I J J1 J2 J3
        
        % create annotation bounding boxes
        [numann,ctlist]=save_bounding_boxes(I0,outpth,nm,numclass);

        numann0=cat(1,numann0,numann);
        ctlist0=cat(1,ctlist0,ctlist);
        
        disp(['  Finished image in ',num2str(round(toc)),' seconds.'])
    end

end