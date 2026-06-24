%% Visualize the whole figure

load('\\Lucie Dequiedt\Kidney Project\Volumes\E11_cleaned.mat')
viewer = viewer3d;
viewer.LightPositionMode ='headlight';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.85;
viewer.CameraZoom = 1.1;
viewer.Denoising = "on";

alpha = ones(size(cleaned))*0.75;
color3 = [97 25 63]/255;
color2 = [255 185 151]/255;
color1 = [246 126 125]/255;

bldcl = imdilate(cleaned==5,strel('sphere',3));
metcl = imdilate(cleaned==1,strel('sphere',3));
mescl = imdilate(cleaned==4,strel('sphere',3));
obj1 = volshow(bldcl, ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = alpha.*bldcl, ...
    Colormap=color1, ... 
    GradientOpacityValue= 0.6);

obj2 = volshow(metcl, ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = alpha.*metcl, ...
    Colormap=color2);

obj3 = volshow(mescl, ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = alpha.*mescl, ...
    Colormap=color3);

