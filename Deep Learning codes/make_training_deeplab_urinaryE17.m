path(path,'base');
warning ('off','all');

% location of annotations and jpg files
pth='\\Lucie Dequiedt\Kidney Project\Mice\E17 mouse\Annotations master model\';
% location of 8-bit tif images corresponding to annotated files
pthim='\\Lucie Dequiedt\E17 mouse\10x\'; 
% location of annotations for model testing
pthtest=[pth,'Testing\'];

% date of model training
nm='03_27_2024_master';

% resolution of tif images (micron / pixel)
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x

% define actions to take per annotation class
% rest kidney WS
WS{1}=[0 0 1];      % remove whitespace if 0, keep only whitespace if 1, keep both if 2
WS{2}=[3];      % add removed whitespace to this class
WS{3}=[1 2 3];      % rename classes accoring to this order 
WS{4}=[1 2 3];      % reverse priority of classes
WS{5}=[];      % delete classes

% color legend for each tissue class
cmap=[245, 223, 171; ... 
    242, 65, 65; ... 
    255 255 255];     

% name of each tissue class
classNames = ["tissue" "kidney" "WS"];

% make check annotation bounding boxes of these tissues for validation
classcheck=[];

% DON'T CHANGE THESE NUMBERS
% number of training tiles
ntrain=10; 
% number of validation tiles
nvalidate=ceil(ntrain/5);
% size of tiles for model training (in pixels)
sxy=1000;

% define path to model based on the model name (don't change this)
pthDL=[pth,nm,'\'];

%% 1 save model metadata
save_model_metadata(pthDL,pthim,WS,nm,umpix,cmap,sxy,classNames,ntrain,nvalidate)

%% 2 load and format annotations from each annotated image
[ctlist0,numann0]=load_annotation_data(pthDL,pth,pthim,classcheck);

% make training & validation tiles for model
create_training_tiles(pthDL,numann0,ctlist0)

% train model
train_segmentation_model(pthDL);

%% test model
pthtestim = '\\Lucie Dequiedt\E17 mouse\Annotations master model\Testing\10x\';

test_segmentation_model(pthDL,pthtest,pthtestim);

%% classify using trained model
pthim = '\\Lucie Dequiedt\Kidney Project\Mice\E17 mouse\10x_reg_chop_K1\';
classify_images(pthim,pthDL);
%%
pthim = '\\Lucie Dequiedt\Kidney Project\Mice\E17 mouse\10x_reg_chop_K2\';
classify_images(pthim,pthDL);

