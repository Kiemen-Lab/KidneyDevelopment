path(path,'base');
warning ('off','all');

% location of annotations and jpg files
pth='\\Lucie Dequiedt\E17 mouse\Annotations kidney model\';
% location of 8-bit tif images corresponding to annotated files
pthim='\\Lucie Dequiedt\E17 mouse\10x\'; 
% location of annotations for model testing
pthtest=['\\Lucie Dequiedt\E17 mouse\Annotations kidney model\Testing\'];

% date of model training
nm='06_20_2025';

% resolution of tif images (micron / pixel)
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x

% define actions to take per annotation class
WS{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 0 0 0 2];      % remove whitespace if 0, keep only whitespace if 1, keep both if 2
WS{2}=[18];      % add removed whitespace to this class
WS{3}=[1 2 3 4 5 6 5 7 8 9 10 11 12 9 1 13 9 14 14 14 15 16 3];      % rename classes accoring to this order 
WS{4}=[10 17 14 12 20 4 6 7 8 9 11 13 15 2 16 19 21 22 23 5 1 3 18];      % reverse priority of classes
WS{5}=[];      % delete classes

% color legend for each tissue class
cmap = [235   134   181; ...  % 1 ureter
    43    76   207; ...       % 2 vein
   134   166   235; ...       % 3 glomerular tuft
    41    70   133; ...       % 4 bowman's capsule
   191   188   109; ...       % 5 distal tubule
   135   214   193; ...       % 6 proximal tubule
   113   191   109; ...       % 7 henle Loop
   235   186   134; ...       % 8 collecting duct
   246   232   250; ...       % 9 loose stroma 
   168   134   235; ...       % 10 developing corpuscle
   144    55   148; ...       % 11 undifferentiated blastema cells 
   100    66   168; ...       % 12 developing nephron 
   145    29    29; ...       % 13 arteries
   255   255   255; ...       % 14 noise 
   235   154   134; ...       % 15 capillary ducts
    64     3     3];          % 16 arterioles
  

% name of each tissue class
classNames = ["urothelium" "vein" "glom_tuft" "bowman_cap"  "dist_tub" "proximal_tub" ...
            "henle_loop" "collecting_duct" "stroma" "dev_corpuscle" "blastema" ...
            "dev_nephron" "arteries" "noise" ... %"capillaries" "noise" ... %"smooth_muscle"...
            "med_coll_duct" "arterioles" "black"];

% make check annotation bounding boxes of these tissues for validation
% classcheck=[1 3 9 10 13 16];
classcheck=[];
% DON'T CHANGE THESE NUMBERS
% number of training tiles
ntrain=12; 
% number of validation tiles
nvalidate=ceil(ntrain/5);
% size of tiles for model training (in pixels)
sxy=1000;

% define path to model based on the model name (don't change this)
pthDL=[pth,nm,'\'];

%% 1 save model metadata
save_model_metadata(pthDL,pthim,WS,nm,umpix,cmap,sxy,classNames,ntrain,nvalidate)

%2 load and format annotations from each annotated image
[ctlist0,numann0]=load_annotation_data(pthDL,pth,pthim,classcheck);
%
% make training & validation tiles for model
create_training_tiles(pthDL,numann0,ctlist0)

% train model
train_segmentation_model(pthDL); 

%% test model
pthtestim = '\\Lucie Dequiedt\E17 mouse\Annotations kidney model\Testing\10x\';
test_segmentation_model(pthDL,pthtest,pthtestim);

%% classify using trained model
pthim = '\\Lucie Dequiedt\Visium kidney\H&E closest slides\E17\10x\';
classify_images(pthim,pthDL);

