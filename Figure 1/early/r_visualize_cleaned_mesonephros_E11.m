load('\\Lucie Dequiedt\Kidney Project\Volumes\E11_cleaned_meso.mat')
tub = cleaned_meso==1;
tub = bwareaopen(tub,200);
tub=imdilate(tub,strel('sphere',3));
color1 = [97 25 63]/255;
alpha = ones(size(tub))*1;
viewer = viewer3d;
viewer.LightPositionMode ='camera-above';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.7; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 1.5;
viewer.Denoising = "on";
obj1 = volshow(tub, ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = alpha.*tub, ...
    Colormap=color1, ... 
    GradientOpacityValue= 0.6);
