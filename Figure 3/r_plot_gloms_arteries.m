%samp = 'Mac_c';
samp = 'E17_K2';
%samp = 'Hum_K3_top';
load(['\\Lucie Dequiedt\Kidney Project\Volumes\',samp,'.mat'])
load(['\\Lucie Dequiedt\Kidney Project\Data for paper\Codes\Figure 4\Glom Volumes\volgloms\',samp,'_labelled_glom.mat']);
mask = imread(['\\Lucie Dequiedt\Kidney Project\Volumes\distK_update\masks\',samp,'.tif']);

art = double(volTA==13);
art = art.*mask;
art = bwareaopen(art,10);
art = imclose(art,strel('sphere',3));
art = imopen(art,strel('sphere',1));
art = bwareaopen(art,50);

arte = double(volTA==16);
arte = arte.*mask;
arte = bwareaopen(arte,10);
arte = imclose(arte,strel('sphere',3));
art = imopen(art,strel('sphere',1));
arte = bwareaopen(arte,50);

%%
viewer =viewer3d(BackgroundColor="white",BackgroundGradient="off",Lighting="on");

volview2 = volshow(volglom~=0,Parent=viewer, ...
    RenderingStyle="VolumeRendering", ...
    Colormap=[134   166   235]/255, ...
        Interpolation = 'nearest', ...
        Alphamap="linear");
volview2 = volshow(arte,Parent=viewer, ...
    RenderingStyle="VolumeRendering", ...
    Colormap=[64     3     3]/255, ...
        Interpolation = 'nearest', ...
        Alphamap="linear");
volview2 = volshow(art,Parent=viewer, ...
    RenderingStyle="VolumeRendering", ...
    Colormap=[145    29    29]/255, ...
        Interpolation = 'nearest', ...
        Alphamap="linear");

