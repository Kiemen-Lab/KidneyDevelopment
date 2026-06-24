load('\\Lucie Dequiedt\Kidney Project\Volumes\Mac_c.mat')
%volTA = volTA(1:2048,1:2048,:);volshow(volTA~=14);
clean = zeros(size(volTA));
%%
stro = volTA==9;
stro = imclose(stro,strel('sphere',2));
stro = imopen(stro,strel('sphere',1));
stro = bwareaopen(stro,250);
clean(stro==1) = 9;


ur1 = volTA==1;
ur = imclose(ur1,strel('sphere',2));
ur = imopen(ur,strel('sphere',1));
ur = bwareaopen(ur,100000);

volTA(ur1&~ur) = 8;
clean(ur) = 1;
%%
glom = volTA==3;
glom = bwareaopen(glom,50);
glom = imclose(glom,strel('sphere',2));
glom = imopen(glom,strel('sphere',2));
glom = bwareaopen(glom,250);
%volshow(glom);
clean(glom) = 3;
%%
dist = volTA==5;
dist = bwareaopen(dist,20);
dist = imclose(dist,strel('sphere',2));
dist = imopen(dist,strel('sphere',1));
dist = bwareaopen(dist,70);
%volshow(dist);
clean(dist) = 5;
%%
prox = volTA==6;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,70);
%volshow(prox);
clean(prox) = 6; 
%%
prox = volTA==7;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,70);
%volshow(prox);
clean(prox) = 7; 
%%
prox = volTA==8;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,1000);
%volshow(prox);
clean(prox) = 8; 
%%
glom = volTA==10;
glom = bwareaopen(glom,50);
glom = imclose(glom,strel('sphere',2));
glom = imopen(glom,strel('sphere',2));
glom = bwareaopen(glom,250);
%volshow(glom);
clean(glom) = 10;
%%
bla = imclose(volTA==11,strel('sphere',3));
bla = bwareaopen(bla,1000000);
clean(bla) = 11;
%volshow(bla);
%%
prox = volTA==12;
prox = bwareaopen(prox,20);
prox = imclose(prox,strel('sphere',2));
prox = imopen(prox,strel('sphere',1));
prox = bwareaopen(prox,70);
%volshow(prox);
clean(prox) = 12;
%%
art = volTA==13;
%volshow(art);
art = bwareaopen(art,10);
art = imclose(art,strel('sphere',2));
%art = imopen(art,strel('sphere',1));
art = bwareaopen(art,70);
%volshow(art);
clean(art) = 13;
%%
med = volTA==15;
med = bwareaopen(med,20);
med = imclose(med,strel('sphere',2));
med = imopen(med,strel('sphere',1));
med = bwareaopen(med,1000);
% volshow(prox);
clean(med==1) = 15;
clean(volTA==16) =16;
clean(volTA==4) =4;
%% Save the cleaned matrix
save('\\Lucie Dequiedt\Kidney Project\Volumes\cleaned\Mac_c.mat','clean');