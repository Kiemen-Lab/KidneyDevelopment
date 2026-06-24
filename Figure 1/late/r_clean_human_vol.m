load('\\Lucie Dequiedt\Kidney Project\Volumes\Hum_K3_top.mat')
%volTA = volTA(1:2048,1:2048,:);volshow(volTA~=14);
cleaned = zeros(size(volTA));
%%
stro = volTA==9;
stro = imclose(stro,strel('sphere',2));
stro = imopen(stro,strel('sphere',1));
stro = bwareaopen(stro,250);
cleaned(stro) = 9;

ur1 = volTA==1;
ur = imclose(ur1,strel('sphere',2));
ur = imopen(ur,strel('sphere',1));
ur = bwareaopen(ur,100000);

volTA(ur1&~ur) = 8;
cleaned(ur) = 1;
%%
glom = volTA==3|volTA==4;
glom = bwareaopen(glom,50);
glom = imclose(glom,strel('sphere',2));
glom = imopen(glom,strel('sphere',2));
glom = bwareaopen(glom,250);
%volshow(glom);
cleaned(glom) = 3;
%%
dist = volTA==5;
dist = bwareaopen(dist,20);
dist = imclose(dist,strel('sphere',2));
dist = imopen(dist,strel('sphere',1));
dist = bwareaopen(dist,70);
%volshow(dist);
cleaned(dist) = 5;
%%
prox = volTA==6;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,70);
%volshow(prox);
cleaned(prox) = 6; 
%%
prox = volTA==7;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,70);
%volshow(prox);
cleaned(prox) = 7; 
%%
prox = volTA==8;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,1000);
%volshow(prox);
cleaned(prox) = 8; 
%%
glom = volTA==10;
glom = bwareaopen(glom,50);
glom = imclose(glom,strel('sphere',2));
glom = imopen(glom,strel('sphere',2));
glom = bwareaopen(glom,250);
%volshow(glom);
cleaned(glom) = 10;
%%
bla = imclose(volTA==11,strel('sphere',3));
bla = bwareaopen(bla,1000000);
cleaned(bla) = 11;
%volshow(bla);
%%
prox = volTA==12;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,70);
%volshow(prox);
cleaned(prox) = 12;
%%
art = volTA==13;
%volshow(art);
art = bwareaopen(art,10);
art = imclose(art,strel('sphere',2));
%art = imopen(art,strel('sphere',1));
art = bwareaopen(art,70);
%volshow(art);
cleaned(art) = 13;
%% clean up medullary collecting duct
med = volTA==15;
med = bwareaopen(med,10000);
med = imclose(med,strel('sphere',2));
med = imfill(med,'holes');
med = imdilate(med,strel('sphere',1));
med = imopen(med,strel('disk',2));
cleaned(med==1) = 15;
cleaned(volTA==16) =16;
cleaned(volTA==4) =4;
%% Save the cleaned matrix
save('\\Lucie Dequiedt\Kidney Project\Volumes\cleaned\Hum_K3_top.mat','cleaned');