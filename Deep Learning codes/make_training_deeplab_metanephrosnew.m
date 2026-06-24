path(path,'base');
warning ('off','all');

% location of annotations and jpg files
pth='\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\Annotate Metanephros\';
% location of 8-bit tif images corresponding to annotated files
pthim='\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\10x_python\'; 
% location of annotations for model testing
pthtest=['\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\Annotate Metanephros\Testing\'];

% date of model training
nm='09_18_2025_metanephrosmouse';

% resolution of tif images (micron / pixel)
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x

% define actions to take per annotation class
WS{1}= [0              0              0                     2            2]; % remove whitespace if 0, keep only whitespace if 1, keep both if 2
WS{2}= 5;                                     % add removed whitespace to nontissue class
WS{3}= [1 2 3 4 5];  % rename classes accoring to this order - cmap order
WS{4}= [4 3 2 1 5];     % reverse priority of classes (bottom to top)
WS{5}= [];                                     % delete classes

% color legend for each tissue class
cmap=[187 219 171;...   % 1  ureteric_tree
    249 239 155;...     % 2  renal_vesicle
    186 143 219;...     % 3  condensed_mesenchyme
    239 181 229;...    % 4  loose_stoma
    255 255 255];    % 5 noise


% name of each tissue class
classNames = ["ureteric_tree" "renal_vesicle" "condensed_mesenchyme" "loose_stoma" "nontissue"];

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
%%
% make training & validation tiles for model
create_training_tiles(pthDL,numann0,ctlist0,0,2,WS{2})

% train model
train_segmentation_model(pthDL);

%% test model
pthtestim = '\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\Annotate Metanephros\Testing\10x\';
test_segmentation_model(pthDL,pthtest,pthtestim);

%% classify using trained model
pthim = '\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\10x_reg_chop\';
classify_images(pthim,pthDL);