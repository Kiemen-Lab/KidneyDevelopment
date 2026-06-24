path(path,'base');
warning ('off','all');

% location of annotations and jpg files
pth='\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\Annotate mesonephros\';
% location of 8-bit tif images corresponding to annotated files
pthim='\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\10x_python\'; 
% location of annotations for model testing
pthtest=['\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\Annotate mesonephros\Testing\'];

% date of model training
nm='10_23_2025_mesonephrosmouse';

% resolution of tif images (micron / pixel)
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x

% define actions to take per annotation class
WS{1}= [2           0                   0               0                  2      0                     2       0]; % remove whitespace if 0, keep only whitespace if 1, keep both if 2
WS{2}= 7;                                     % add removed whitespace to nontissue class
WS{3}= [1 2 3 4 5 6 7 8];  % rename classes accoring to this order - cmap order
WS{4}= [5 3 2 4 6 1 8 7];     % reverse priority of classes (bottom to top)
WS{5}= [1 2 3 6 8];      % delete classes

% color legend for each tissue class
cmap=[157 202 194;...    % 4  mesonephric_duct
    224 187 228; ...     % 5  ECM
    255 255 255];        % nontissue


% name of each tissue class
classNames = ["Mesonephric_Duct" "ECM" "nontissue"];

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
create_training_tiles(pthDL,numann0,ctlist0,0,3,3)

% train model
train_segmentation_model(pthDL);

%% test model
pthtestim = '\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\Annotate mesonephros\Testing\10x\';
test_segmentation_model(pthDL,pthtest,pthtestim);

%% classify using trained model
pthim = '\\Lucie Dequiedt\Kidney Project\Mice\E11 mouse\10x_reg_chop\';
classify_images(pthim,pthDL);