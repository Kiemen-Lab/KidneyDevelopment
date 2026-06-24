warning ('off','all');
pth='\lucie\Fetal Kidney\kidney\annotate\';
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x
pthim='\collaborations\monkey\kidney\whole sample\10x\';
pthclassify_a=[pthim,'chopped\a\']; % path to tif images to classify chop a
pthclassify_b=[pthim,'chopped\b\']; % path to tif images to classify chop b
pthclassify_c=[pthim,'chopped\c\']; % path to tif images to classify chop c
pthclassify_d=[pthim,'chopped\d\']; % path to tif images to classify chop d

pthTA=[pthim,'TA\'];
nm='deeplab_07_18_glom_comb_WS_NOTremoved_BCcombinetile_AP';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIRST MODEL
% 1 ureter
% 2 vein
% 3 glomerulus 
% 4 Bowman's capsule 
% 5 macula densa
% 6 proximal convoluted tubule
% 7 distal convoluted tubule
% 8 Henle loop
% 9 collecting duct
% 10 loose stroma 
% 11 developing corpuscle
% 12 undifferentiated blastema cells
% 13 developing nephron
% 14 conjunctive capsule 
% 15 renal pelvis
% 16 arteries
% 17 capillaries
% 18 noise 
% 19 nerves 
%%
%define actions to take per annotation class
% WS{1}=[0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2];               % remove whitespace if 0, keep only whitespace if 1, keep both if 2
% WS{2}=18;                                                     % add removed whitespace to this class - NOISE
% % combine 7/8/10 - whitespace;  combine 2/12 - duct & metaplasia
% WS{3}=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19];      % rename classes accoring to this order 
% WS{4}=[10 14 12 4 6 7 8 9 11 13 15 2 16 19 5 1 3 17 18];      % reverse priority of classes
% WS{5}=[];                                                     % delete classes
% numclass=length(unique(WS{3}));
% sxy=750;
% pthDL=[pth,nm,'\'];
% nblack=numclass+1;
%%

% cmap=[108 75 189;... % 1 ureter
%     43 76 207;...    % 2 vein
%     96 166 5;...     % 3 glomerulus
%     215 222 91;...   % 4 Bowman's capsule
%     149 35  184;...  % 5 macula densa
%     222 187 91;...   % 6 proximal convoluted tubule
%     222 91 91;...    % 7 distal convoluted tubule
%     222 135 91;...   % 8 Henle loop
%     189 75 162;...   % 9 collecting duct
%     235 201 245; ... % 10 loose stroma
%     110 129 240; ... % 11 developing corpuscle
%     209 110 240; ... % 12 undifferentiated blastema cells
%     130 110 240; ... % 13 developing nephron
%     145 4 184; ...   % 14 conjunctive capsule
%     157 75 189; ...  % 15 renal pelvis
%     237 36 36; ...   % 16 arteries
%     122 79 143; ...  % 17 capillaires
%     255 255 255; ... % 18 noise
%     104 86 110;];    % 19 nerves
% 
% classNames=["ureter" "vein" "glomerulus" "bowman_cap" "macula_densa" "proximal_tub" ...
%             "distal_tub" "henle_loop" "collecting_duct" "stroma" "dev_corpuscle" "blastema" ...
%             "dev_nephron" "capsule" "pelvis" "arteries" "capillaries" "noise" "nerves" "black"];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW MODEL 
% 1 ureter                              1
% 2 vein                                2 
% 3 glomerulus                          3  
% 4 Bowman's capsule                    4
% 5 macula densa                        5
% 6 proximal convoluted tubule          6
% 7 distal convoluted tubule            5
% 8 Henle loop                          7
% 9 collecting duct                     8
% 10 loose stroma                       9
% 11 developing corpuscle               10  
% 12 undifferentiated blastema cells    11 
% 13 developing nephron                 12
% 14 conjunctive capsule                13
% 15 renal pelvis                       1
% 16 arteries                           14
% 17 capillaries                        9
% 18 noise                              15
% 19 nerves                             13
% 20 smooth muscle                      13
% 21 capillary duct                     16
% 22 arterioles                         17   
% 23 fetal glomerular tuft              18
%%
% define actions to take per annotation class
WS{1}=[0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 0 0 0 2];               % remove whitespace if 0, keep only whitespace if 1, keep both if 2
WS{2}=18;                                                    % add removed whitespace to this class - NOISE
% % combine 7/8/10 - whitespace;  combine 2/12 - duct & metaplasia
WS{3}=[1 2 3 4 5 6 5 7 8 9 10 11 12 13 1 14 9 15 13 13 16 17 3];      % rename classes accoring to this order 
WS{4}=[10 17 14 12 20 4 6 7 8 9 11 13 15 2 16 19 21 22 23 5 1 3 18];      % reverse priority of classes
WS{5}=[];                                                     % delete classes
numclass=length(unique(WS{3}));
sxy=750;
pthDL=[pth,nm,'\'];
nblack=numclass+1;
nwhite=15;
cmap=[235 134 181; ...    % 1 ureter
    43 76 207;...         % 2 vein
    134 166 235;...       % 3 glomerulus  
    41 70 133;...         % 4 Bowman's capsule 
    %149 35  184;...      % 5 macula densa
    135 214 193; ...      % 6 proximal convoluted tubule 
    191 188 109; ...      % 7 distal convoluted tubule 5 107 166
    113 191 109; ...      % 8 Henle loop 
    235 186 134; ...      % 9 collecting duct
    246 232 250; ...      % 10 loose stroma
    168 134 235; ...      % 11 developing corpuscle
    144 55 148; ...       % 12 undifferentiated blastema cells 
    100 66 168; ...       % 13 developing nephron 
    225 159 245; ...      % 14 conjunctive capsule
    %235 154 134; ...     % 15 renal pelvis
    145 29 29; ...        % 16 arteries
    45 5 56; ...          % 17 capillaires
    255 255 255; ...      % 18 noise
    %104 86 110; ...      % 19 nerves
    191 159 174];         % 20 smooth muscle


%classNames=["ureter" "vein" "glomerulus" "bowman_cap" "macula_densa" "proximal_tub" ...
         %    "distal_tub" "henle_loop" "collecting_duct" "stroma" "dev_corpuscle" "blastema" ...
         %    "dev_nephron" "capsule" "pelvis" "arteries" "capillaries" "noise" "nerves" "smooth_muscle" "black"];

% classNames=["pelvis" "vein" "glomerulus" "bowman_cap"  "proximal_tub" ...
%             "distal_tub" "henle_loop" "collecting_duct" "stroma" "dev_corpuscle" "blastema" ...
%             "dev_nephron" "capsule" "arteries" "capillaries" "noise" "smooth_muscle" "black"];
%         
classNames=["urothelium" "vein" "glom_tuft" "bowman_cap"  "dist_tub" "proximal_tub" ...
            "henle_loop" "collecting_duct" "stroma" "dev_corpuscle" "blastema" ...
            "dev_nephron" "conj_cap" "arteries" "noise" ... %"capillaries" "noise" ... %"smooth_muscle"...
            "cap_duct" "arterioles" "black"]; %"fet_glom_tuft" "black" ];

%% load and format annotations for each image
% for each annotation file
imlist=dir([pth,'*xml']);
numann0=[];ctlist=[];
for kk=1:length(imlist)
    % set up names
    imnm=imlist(kk).name(1:end-4);tic;
    disp(['Image ',num2str(kk),' of ',num2str(length(imlist)),': ',imnm])
    outpth=[pth,'data\',imnm,'\'];
    if ~exist(outpth,'dir');mkdir(outpth);end
    matfile=[outpth,'annotations.mat'];
    
    % skip if file hasn't been updated since last load
    dm='';bb=0;date_modified=imlist(kk).date;
    if exist(matfile,'file');load(matfile,'dm','bb');end
    if contains(dm,date_modified) && bb==1
        disp('  annotation data previously loaded')
        load([outpth,'annotations.mat'],'numann','ctlist0');
        numann0=[numann0;numann];ctlist=[ctlist;ctlist0];
        continue;
    end
    
    % 1 read xml annotation files and saves as mat files
    load_xml_file(outpth,[pth,imnm,'.xml'],date_modified);
    
     % 2 fill annotation outlines and delete unwanted pixels
    [I,TA,pthTA]=calculate_tissue_space_082(pthim,imnm);
    J=fill_annotations_file(I,outpth,WS,umpix,TA); 
    J=edit_individual_classes_082(J,[outpth,'view_annotations.tif']);
    figure(9);
        subplot(1,2,1);imshow(uint8(I))
        subplot(1,2,2);imagesc(J);axis equal;axis off
        ha2=get(gcf,'children');linkaxes(ha2);
    
    
    % create original and color corrected annotation bounding boxes
    [numann,ctlist0]=annotation_bounding_boxes(I,outpth,nm,numclass);
    numann0=[numann0;numann];ctlist=[ctlist;ctlist0];
    toc;
end


%% combine tiles 
%path('\\babyserverdw3\ashley kiemen\PanIN Modelling Package\Fetal Kidney\classification\')
numann=numann0;
% make training tiles
ty='training\';obg=[pthDL,ty,'big_tiles\'];
while length(dir([obg,'*tif']))<14
    numann=combine_tiles_big(numann0,numann,ctlist,nblack,pthDL,ty,sxy,1);
    disp([num2str(length(dir([obg,'*tif']))/2),' images complete'])
end
% make validation tiles
ty='validation\';obg=[pthDL,ty,'big_tiles\'];
while length(dir([obg,'*tif']))<2
    numann=combine_tiles_big(numann0,numann,ctlist,nblack,pthDL,ty,sxy,1);
    disp([num2str(length(dir([obg,'*tif']))/2),' images complete'])
end

%% build model
train_deeplab(pthDL,1:numclass+1,sxy,classNames);


 %% classify all images

% % pth=[pthclassify_a,pthclassify_b,pthclassify_c,pthclassify_d]
% % deeplab_classification(pthclassify_a,pthDL,sxy,nm,cmap,nblack,16);
% 
deeplab_classification(pthclassify_a,pthDL,sxy,nm,cmap,nblack,nwhite);
% deeplab_classification(pthclassify_b,pthDL,sxy,nm,cmap,nblack,nwhite);
% deeplab_classification(pthclassify_c,pthDL,sxy,nm,cmap,nblack,nwhite);
% deeplab_classification(pthclassify_d,pthDL,sxy,nm,cmap,nblack,nwhite);
% 
%  % register classified images
% pthdata_a='\\motherserverdw\Expansion\ashleyserverdw\collaborations\monkey\kidney\whole sample\1x\a\registered\elastic registration\save_warps\';
% scale=8;
% pthim_a=[pthclassify_a,'classification_',nm,'\'];
% save_images_elastic2(pthim_a,pthdata_a,scale,nwhite);
% 
% pthdata_b='\\motherserverdw\Expansion\ashleyserverdw\collaborations\monkey\kidney\whole sample\1x\b\registered\elastic registration\save_warps\';
% scale=8;
% pthim_b=[pthclassify_b,'classification_',nm,'\'];
% save_images_elastic2(pthim_b,pthdata_b,scale,nwhite);
% 
% pthdata_c='\\motherserverdw\Expansion\ashleyserverdw\collaborations\monkey\kidney\whole sample\1x\c\registered\elastic registration\save_warps\';
% scale=8;
% pthim_c=[pthclassify_c,'classification_',nm,'\'];
% save_images_elastic2(pthim_c,pthdata_c,scale,nwhite);
% 
% pthdata_d='\\motherserverdw\Expansion\ashleyserverdw\collaborations\monkey\kidney\whole sample\1x\d\registered\elastic registration\save_warps\';
% scale=8;
% pthim_d=[pthclassify_d,'classification_',nm,'\'];
% save_images_elastic2(pthim_d,pthdata_d,scale,nwhite);
