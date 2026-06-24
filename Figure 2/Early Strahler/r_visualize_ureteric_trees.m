load('\\Lucie Dequiedt\Kidney Project\Volumes\E11_cleaned_meta.mat')
load("\\Lucie Dequiedt\Kidney Project\Volumes\metanephros\E11_meta.mat")
skeleton = information.skeleton{1};
strahlerLabels = information.generationLabels{1};
vol = cleaned_meta;
clearvars cleaned_meta


cmap = [1 1 1; ... % background
    0.0196    0.5412    1.0000; ... % blue
    0.0196    0.9882    0.5373; ...  % cyan
   1 1 0; ...  % green
 1 0.5 0; ... % orange
1 0 0];    % red

%% create tree object
kidneys = vol~=0;
kidneys = bwareaopen(kidneys,10000);
kidneys = bwlabeln(kidneys);
caps = double(vol==3).*kidneys;
inside = imdilate(caps,strel('sphere',5));
inside = bwareaopen(inside,100000);
inside = imdilate(inside,strel('sphere',10));
% Load or create your binary volume
tree = double(vol==1).*kidneys;
tree = tree.*double(inside);

tree = tree==1;
tree = bwareaopen(tree,100);
tree = imclose(tree,strel('sphere',5));
tree = imfill(tree,'holes');
t = zeros(size(tree));
for i=1:size(tree,3)
    t(:,:,i) = imfill(tree(:,:,i),'holes');
end
tree = t;
%% visualize tree only
viewer = viewer3d;
viewer.LightPositionMode ='camera-above';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 2.5;
viewer.Denoising = "on";

alpha = ones(size(tree(:,:,:)));
color1 = [255 185 151]/255;
%
% obj1 = volshow(tree, ...
%     Parent=viewer, ...
%     RenderingStyle="GradientOpacity", ...
%     LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
%     SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
%     AlphaData = alpha.*tree, ...
%     Colormap=color1, ... 
%     GradientOpacityValue= 0.6);

% visualize tree + skeleton 
viewer = viewer3d;
viewer.LightPositionMode ='camera-above';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 2.5;
viewer.Denoising = "on";

alpha = ones(size(tree(:,:,:)));
color1 = [184, 178, 178]/255;
volshow(tree,Parent=viewer, ...
    RenderingStyle="LightScattering", ...
    Colormap=color1, ...
    Alphamap="linear");

skeletonLabeled = zeros(size(skeleton));
skeletonLabeled(skeleton) = strahlerLabels;
skl = imdilate(skeletonLabeled,strel('sphere',2));


volshow(skl,Parent=viewer, ...
    Interpolation='nearest', ...
    RenderingStyle="LightScattering", ...
    Colormap=cmap(1:max(skl(:)+1),:), ...
    Alphamap=[0 ones(1,5)]);



%%
clearvars 
load('\\Lucie Dequiedt\Kidney Project\Volumes\MA40_cleaned_meta.mat')
load("\\Lucie Dequiedt\Kidney Project\Volumes\metanephros\MA40_meta.mat")
skeleton = information.skeleton{1};
strahlerLabels = information.generationLabels{1};

vol = cleaned_meta;
clearvars cleaned_meta
%% create tree object
kidneys = vol~=0;
kidneys = bwareaopen(kidneys,10000);
kidneys = bwlabeln(kidneys);
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
% Load or create your binary volume
tree = double(vol==1).*kidneys;
tree = tree.*double(inside);

tree = tree==1;
tree = bwareaopen(tree,100);
tree = imfill(tree,'holes');
tree = imclose(tree,strel('sphere',2));tree = imopen(tree,strel('sphere',2));
tree = bwareaopen(tree,10000);
t = zeros(size(tree));
for i=1:size(tree,3)
    t(:,:,i) = imfill(tree(:,:,i),'holes');
end
tree = t;
clearvars t 
%% visualize tree only
cmap = [1 1 1; ... % background
    0.0196    0.5412    1.0000; ... % blue
    0.0196    0.9882    0.5373; ...  % cyan
   1 1 0; ...  % green
 1 0.5 0; ... % orange
1 0 0];    % red
viewer = viewer3d;
viewer.LightPositionMode ='camera-above';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 2.5;
viewer.Denoising = "on";

alpha = ones(size(tree(:,:,:)));
color1 = [255 185 151]/255;
% obj1 = volshow(tree, ...
%     Parent=viewer, ...
%     RenderingStyle="GradientOpacity", ...
%     LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
%     SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
%     AlphaData = alpha.*tree, ...
%     Colormap=color1, ... 
%     GradientOpacityValue= 0.6);

% visualize tree + skeleton 
viewer = viewer3d;
viewer.LightPositionMode ='camera-above';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 2.5;
viewer.Denoising = "on";

alpha = ones(size(tree(:,:,:)));
color1 = [184, 178, 178]/255;
volshow(tree,Parent=viewer, ...
    RenderingStyle="LightScattering", ...
    Colormap=color1, ...
    Alphamap=[0 0.6]);

skeletonLabeled = zeros(size(skeleton));
skeletonLabeled(skeleton) = strahlerLabels;
skl = imdilate(skeletonLabeled,strel('sphere',2));


volshow(skl,Parent=viewer, ...
    Interpolation='nearest', ...
    RenderingStyle="LightScattering", ...
    Colormap=cmap(1:max(skl(:)+1),:), ...
    Alphamap=[0 ones(1,5)]);

%%

% Define your key colors
cmap_key = [0.0196    0.5412    1.0000; ... % blue
    0.0196    0.9882    0.5373; ...  % cyan
   1 1 0; ...  % yellow
 1 0.5 0; ... % orange
1 0 0];    % red

% Create positions for each key color (1 to 5 for Strahler orders)
key_positions = 1:5;  % Changed from 0:5 to match 5 colors

% Create smooth colormap with 256 levels
n_colors = 256;
query_positions = linspace(1, 5, n_colors);  % Changed from 0 to 5
cmap_smooth = interp1(key_positions, cmap_key, query_positions, 'linear');

% Create a figure just for the colorbar
figure('Position', [100, 100, 200, 400]);
imagesc([1:5]');  % Vertical array (with transpose)
colormap(cmap_smooth);
c = colorbar;  % Vertical colorbar (default)
c.Ticks = 1:5;  % Changed from 0:5
c.TickLabels = {'1', '2', '3', '4', '5'};
c.Label.String = 'Strahler Order';
c.Label.FontSize = 12;
axis off;
caxis([1 5]);  % Set color axis limits to match data