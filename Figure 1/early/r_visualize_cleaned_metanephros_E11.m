load('\\Lucie Dequiedt\Kidney Project\Volumes\E11_cleaned_meta.mat')
%%  ok
tub = cleaned_meta==1;
tub = bwareaopen(tub,1000);
tub = imclose(tub,strel('sphere',5));
tub = imfill(tub,'holes');
t = zeros(size(tub));
for i=1:size(tub,3)
    t(:,:,i) = imfill(tub(:,:,i),'holes');
end
tub = t;
tub = imdilate(tub,strel('sphere',2));
clearvars t 
%%  ok
rn = cleaned_meta==2;
rn = imclose(rn,strel('sphere',2));
rn = imopen(rn,strel('sphere',1));
rn = bwareaopen(rn,200);
rn = imdilate(rn,strel('sphere',2));
%%  ok
cm = cleaned_meta==3;
cm = imclose(cm,strel('sphere',1));
cm = imopen(cm,strel('sphere',1));
cm=bwareaopen(cm,10000);
%cm=imfill(cm,'holes');
%%  whole tubules and renal vesicles
viewer = viewer3d;
viewer.LightPositionMode ='camera-above';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 2.5;
viewer.Denoising = "on";

alpha = ones(size(tub(:,:,80:130)))*0.75;
color1 = [255 185 151]/255;
color2 = [135 76 57]/255;
color3 = [193 174 167]/255;
obj1 = volshow(tub(:,:,:), ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = (ones(size(tub(:,:,:)))*0.75).*tub(:,:,:), ...
    Colormap=color1, ... 
    GradientOpacityValue= 0.6);

obj2 = volshow(rn(:,:,:), ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = (ones(size(tub(:,:,:)))*0.75).*rn(:,:,:), ...
    Colormap=color2);

%% %% half way through 

viewer = viewer3d;
viewer.LightPositionMode ='headlight';
viewer.BackgroundColor = 'white';
viewer.BackgroundGradient = 'off';
viewer.CameraPosition = [306.9720  598.9130 -680.6713];
viewer.LightColor = [255, 250, 214]/255;
viewer.DiffuseLight = 0.65; % 0.7 for big both kidneys, 0.65 for small subset 
viewer.CameraZoom = 2.5;
viewer.Denoising = "on";

color1 = [255 185 151]/255;
color2 = [135 76 57]/255;
color3 = [193 174 167]/255;
stt =65;
endd =140;

obj3 = volshow(cm(:,:,stt:endd), ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = (ones(size(tub(:,:,stt:endd)))*0.75).*cm(:,:,stt:endd), ...
    Colormap=color3);
obj1 = volshow(tub(:,:,stt:endd), ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = (ones(size(tub(:,:,stt:endd)))*0.75).*tub(:,:,stt:endd), ...
    Colormap=color1, ... 
    GradientOpacityValue= 0.6);

obj2 = volshow(rn(:,:,stt:endd), ...
    Parent=viewer, ...
    RenderingStyle="GradientOpacity", ...
    LightScatteringQuality= 0.2, ... % affects visualization only when RenderingStyle is LightScattering, larger value = more realistic rendering ([0,1])
    SpecularReflectance=0.4, ... % This value controls the amount of light reflected by the volume. Increase the reflectance to make the volume appear shinier
    AlphaData = (ones(size(tub(:,:,stt:endd)))*0.75).*rn(:,:,stt:endd), ...
    Colormap=color2);

