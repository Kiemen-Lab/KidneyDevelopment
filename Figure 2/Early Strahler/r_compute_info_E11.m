load('\\Lucie Dequiedt\Kidney Project\Volumes\E11_cleaned_meta.mat')

vol = cleaned_meta;
clearvars cleaned_meta

kidneys = vol~=0;
kidneys = bwareaopen(kidneys,10000);
kidneys = bwlabeln(kidneys);

%% Compute size of tips and total kidney volume
caps = double(vol==3).*kidneys;
inside = imdilate(caps,strel('sphere',5));
inside = bwareaopen(inside,100000);
inside = imdilate(inside,strel('sphere',10));


tips = double(vol==2).*kidneys;
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
tree = tree.*double(inside);
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

[information.numGenerations{1}, information.generationLabels{1},information.skeleton{1}] = compute_Strahler(tree,1);
[information.numGenerations{2}, information.generationLabels{2},information.skeleton{2}] = compute_Strahler(tree,2);

save(['\\Lucie Dequiedt\Kidney Project\Volumes\metanephros\E11_meta.mat'],'information');
