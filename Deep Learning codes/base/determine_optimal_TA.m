function determine_optimal_TA(pthim,numims)
close all;

if nargin==0
    pthim='\\10.17.97.73\kiemen-lab-data\tissue scans\Pancreas Metastases\20231011_Lucie_Pancreas_Liver\Liv-17\10x\';
    numims=5;
end
    
    CT0=205;
    CTA=CT0+5;
    CTC=CT0-5;
    szz=600;

    % choose images to load
    imlist=dir([pthim,'*tif']);
    numims=min([length(imlist) numims]);
    kk=randperm(length(imlist),numims);
    imlist=imlist(kk);
    disp(['Evaluating ',num2str(numims),' randomly selected images to choose a good whitespace detection...'])

    % check if we'vea already done this
    outpth=[pthim,'TA\'];
    if ~isfolder(outpth);mkdir(outpth);end
    if exist([outpth,'TA_cutoff.mat'],'file')
        disp('  Optimal cutoff already chosen, skip this step');
        return;
    end

    cts=zeros([1 length(imlist)]);
    for b=1:length(imlist)
        % load image
        nm=imlist(b).name;
        disp(['  Loading image ',num2str(b),' of ',num2str(numims),': ',nm])
        im0=imread([pthim,nm]);
        disp('    Image loaded')
        
        % crop a region from the image
        do_again=1;
        while do_again==1
            waitfor(msgbox(["An image will appear";" ";...
                "Click on a location at the edge of tissue and whitespace"]));
            h=figure(122);imshow(im0);
            [x,y]=ginput(1);
            try close(h);catch;end
            x=round(x);
            y=round(y);
            xx=intersect(1:size(im0,2),x-szz:x+szz);
            yy=intersect(1:size(im0,1),y-szz:y+szz);
            
            % determine if the cropped region is good or not
            im=im0(yy,xx,:);
            h=figure;imshow(im);
            answer = questdlg('Is this a good location to evaluate tissue and whitespace detection?',...
	            'Evaluate cropped region','Looks good', 'No, select a new location','Looks good');
    
            % Handle response
            switch answer
                case 'Looks good'
                    do_again=0;
                case 'No, select a new location'
                    do_again=1;
            end
            try close(h);catch;end
        end

        found_CT0=0;
        while found_CT0==0
            img=im(:,:,2);title('H&E image')
            imA=img>CTA;imA=bwareaopen(imA,10);imA=~bwareaopen(~imA,10);
            imB=img>CT0;imB=bwareaopen(imB,10);imB=~bwareaopen(~imB,10);
            imC=img>CTC;imC=bwareaopen(imC,10);imC=~bwareaopen(~imC,10);
            
            figure('units','normalized','outerposition',[0 0 1 1]);h=gcf;
                subplot(2,3,2);imshow(im);title('H&E image')
                subplot(2,3,4);imshow(labeloverlay(im,imA,'Color',[0 0 0;1 0.6 0.6]));title('IMAGE A: keep more tissue')
                subplot(2,3,5);imshow(labeloverlay(im,imB,'Color',[0 0 0;1 0.6 0.6]));title('IMAGE B: medium')
                subplot(2,3,6);imshow(labeloverlay(im,imC,'Color',[0 0 0;1 0.6 0.6]));title('IMAGE C: keep more whitespace')
                ha=get(gcf,'children');linkaxes(ha);
    
            waitfor(msgbox(["Zoom in to evaluate the potential whitespace detection images";" ";...
                "Determine which image is best, then click spacebar to continue";" ";...
                "In the images, black pixels are tissue, white pixels are whitespace"]));
        
            zoom on;
            disp('    Press spacebar to return to annotation')
            waitfor(gcf,'CurrentCharacter',' ');
            
            answer = questdlg('Does one of the images look good?','Decide on images',...
	            'Yes', 'No, keep more tissue','No, keep more whitespace','Yes');
            % Handle response
            switch answer
                case 'Yes'
                    found_CT0=1;
                case 'No, keep more tissue'
                    CT0=CT0+10;CTA=CT0+5;CTC=CT0-5;
                    try close(h);catch;end
                case 'No, keep more whitespace'
                    CT0=CT0-10;CTA=CT0+5;CTC=CT0-5;
                    try close(h);catch;end
            end
        end

        answer = questdlg('Which image is best','Choose a whitespace detection cutoff. If all images look good, choose Image B',...
            'IMAGE A', 'IMAGE B','IMAGE C','IMAGE C');
        try close(h);catch;end
        % Handle response
        switch answer
            case 'IMAGE A'
                disp('    chose Image A')
                cts(b)=CTA;
            case 'IMAGE B'
                disp('    chose Image B')
                cts(b)=CT0;
            case 'IMAGE C'
                disp('    chose Image C')
                cts(b)=CTC;
        end
        
        close all;
    end

    save([outpth,'TA_cutoff.mat'],'imlist','cts')
end


