load('\\Lucie Dequiedt\Kidney Project\Volumes\E17_K2.mat')
vol1 = volTA;
clearvars volTA
im=sum(vol1,3);
im=im/max(im(:));
im2=sum(vol1,3);im2=im2/max(im2(:));figure;imagesc(im2);
mask=freeform_annotation(im); % annotate the region you want to KEEP (exclude the noise you want deleted)
vol1=double(vol1).*mask;
%%

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
    64     3     3; ...       % 16 arterioles
    255 255 255];    

viewer = viewer3d;
viewer.LightPositionMode ='auto';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65;
viewer.CameraZoom = 1.1;
viewer.Denoising = "on";
% vol1 = vol(:,1:2048,1:100);
alpha = ones(size(vol1))*0.5;
A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
tform = affinetform3d(A);
%
voll = vol1(:,:,1:75);
for i=[1 2 3 4 5 6 7 8 10 11 12 13 15 16]
    p = logical(vol1==i);
    p=bwareaopen(p,1000);
    p=imdilate(p,strel('sphere',2));
    obj1 = volshow(p(:,:,1:119), ...
        Parent=viewer, ...
        RenderingStyle="GradientOpacity", ...
        LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
        SpecularReflectance=0.2, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
        Colormap=cmap(i,:)/255, ... 
        GradientOpacityValue= 0.6, ...
        Transformation=tform);
      %  AlphaData = alpha.*logical(vol1), ...
end