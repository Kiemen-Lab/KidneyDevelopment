pth = '\\Lucie Dequiedt\Kidney Project\Volumes\tubule_volumes_cleaned\';

samp = 'Mac_c.mat';
load([pth,samp]);

% Tubule colors
cmap = [191   188   109; ...       % 5 distal tubule
        135   214   193; ...       % 6 proximal tubule
        113   191   109; ...       % 7 henle Loopsa
        235   186   134]/255;      % 8 collecting duct

%%
viewer = viewer3d(BackgroundColor="white",BackgroundGradient="off",Lighting="on",AmbientLight=0.7,LightPosition=[1500 0 -1000]);

col = volTA ==8;
col = bwareaopen(col,1000);

volview2 = volshow(col(:,:,:),Parent=viewer, ...
    Colormap=cmap(4,:), ...
    Interpolation = 'bilinear',Alphamap="linear");

henle = volTA==7;
henle = bwareaopen(henle(:,:,:),2000);
volview2 = volshow(henle,Parent=viewer, ...
    Colormap=cmap(3,:), ...
    Interpolation = 'bilinear',Alphamap="linear");

% DISTAL TUBULE - with transparency
volview2 = volshow(volTA(:,:,:)==5,Parent=viewer,...
    Colormap=cmap(1,:), ...
    Interpolation = 'bilinear',Alphamap='linear');  % Adjust transparency (0-1)

% PROXIMAL TUBULE - with transparency
volview2 = volshow(volTA(:,:,:)==6,Parent=viewer, ...
    Colormap=cmap(2,:), ...
    Interpolation = 'bilinear',Alphamap='linear');  % Adjust transparency (0-1)
%%
pth = '\\10.99.134.183\kiemen-lab-data\Lucie Dequiedt\Kidney Project\Volumes\';
load([pth,samp]);

ur = bwareaopen(volTA==1,100000);
ur = imdilate(ur,strel('sphere',7));

mc = bwareaopen(volTA==15,500);
mc = imclose(mc,strel('sphere',3));

%%
viewer = viewer3d(BackgroundColor="white",BackgroundGradient="off",Lighting="on",AmbientLight=0.5,LightPosition=[1500 0 -1000]);

volview2 = volshow(ur(:,:,:),Parent=viewer, ...
    Colormap=[224, 168, 194]/255, ...
    Interpolation = 'bilinear',Alphamap="linear");
volview2 = volshow(mc(:,:,:),Parent=viewer, ...
    Colormap=[235   154   134]/255, ...
    Interpolation = 'bilinear',Alphamap="linear");

volview2 = volshow(col(:,:,:),Parent=viewer, ...
    Colormap=cmap(4,:), ...
    Interpolation = 'bilinear',Alphamap="linear");