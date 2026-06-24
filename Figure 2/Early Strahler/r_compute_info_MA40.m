load('\\Lucie Dequiedt\Kidney Project\Volumes\MA40_cleaned_meta.mat')

vol = cleaned_meta;
clearvars cleaned_meta

% 103 -> 228
kidneys = vol~=0;
kidneys = bwareaopen(kidneys,10000);
kidneys = bwlabeln(kidneys);

%% Compute size of tips and total kidney volume
caps = double(vol==3).*kidneys;


tips = double(vol==2).*kidneys;

inside = kidneys;
inside(:,:,1:102) =0;
inside(:,:,230:325) =0;
for i=365:400
    im = kidneys(:,:,i);
    im = imfill(im,'holes');
    im = bwareafilt(logical(im),1);
    inside(:,:,i) = double(im)*2;
end
inside = bwareaopen(inside,10000);
% inside = imdilate(caps|tips,strel('sphere',5));
% inside = bwareaopen(inside,100000);
% inside = imdilate(inside,strel('sphere',25));

tips = tips.*double(inside);

information = table;
information.Name = ['K1';'K2'];

v = regionprops3(tips==1,'Volume');
information.TipsVolumes{1} = v.Volume;
information.TotVol{1} = sum(sum(sum(inside.*kidneys==1)));
v = regionprops3(tips==2,'Volume');
information.TipsVolumes{2} = v.Volume;
information.TotVol{2} = sum(sum(sum(inside.*kidneys==2)));

%% Volume of ureteric tree
tree = double(vol==1).*kidneys;

information.TreeVol{1} = sum(sum(sum(inside.*tree==1)));
information.TreeVol{2} = sum(sum(sum(inside.*tree==2)));

%% Volume of stroma
stro = double(vol==4).*kidneys;

information.StrVol{1} = sum(sum(sum(inside.*stro==1)));
information.StrVol{2} = sum(sum(sum(inside.*stro==2)));

%% Volume of cap
information.CapVol{1} = sum(sum(sum(inside.*caps==1)));
information.CapVol{2} = sum(sum(sum(inside.*caps==2)));

%% Branching order of tree

% Load or create your binary volume
tree = double(vol==1).*kidneys;
tree = tree.*double(inside);
% tree = double(tree).*double(volTA(:,:,1:545)==2);

clean =1;
[information.numGenerations{1}, information.generationLabels{1},information.skeleton{1}] = compute_Strahler(tree,1,clean);
[information.numGenerations{2}, information.generationLabels{2},information.skeleton{2}] = compute_Strahler(tree,2,clean,[330 281 362]);

save(['\\Lucie Dequiedt\Kidney Project\Volumes\metanephros\MA40_meta.mat'],'information');
