load('\\Lucie Dequiedt\Kidney Project\Volumes\E17_K1.mat')

mask = imread(['\\Lucie Dequiedt\Kidney Project\Volumes\distK_update\masks\E17_K1.tif']);

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


volview2 = volshow(arte,Parent=viewer, ...
    RenderingStyle="CinematicRendering", ...
    Colormap=[64     3     3]/255, ...
        Interpolation = 'nearest', ...
        Alphamap="linear");
volview2 = volshow(art,Parent=viewer, ...
    RenderingStyle="CinematicRendering", ...
    Colormap=[145    29    29]/255, ...
        Interpolation = 'nearest', ...
        Alphamap="linear");
