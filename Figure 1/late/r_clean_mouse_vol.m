load('\\Lucie Dequiedt\Kidney Project\Volumes\E17_K1.mat')
cleaned = zeros(size(volTA));

%%
stro = volTA==9;
stro = imclose(stro,strel('sphere',2));
stro = imopen(stro,strel('sphere',1));
stro = bwareaopen(stro,250);
cleaned(stro==1) = 9;

%% clean up urothelium
ur = volTA==1;
ur = imclose(ur,strel('disk',3));
%ur = imopen(ur,strel('disk',1));
ur = bwareaopen(ur,300000);
cleaned(ur==1) = 1;

%% clean up veins
ves = volTA==2;
ves = bwareaopen(ves,250);
ves = imclose(ves,strel('sphere',2));
ves = imfill(ves,'holes');
ves = bwareaopen(ves,10000);
ves = imdilate(ves,strel('sphere',1));
%ves = imopen(ves,strel('sphere',1));
ves = imopen(ves,strel('disk',1));
cleaned(ves==1) = 2;

%% clean up glomeruli and bowmans
glom = volTA==3;
glom = imclose(glom,strel('sphere',2));
glom = imopen(glom,strel('sphere',1));
glom = bwareaopen(glom,250);
glom = imdilate(glom,strel('disk',2));
glom = imopen(glom,strel('disk',1));
cleaned(glom==1) = 3;

%% clean up distal
dist = volTA==5;
dist = bwareaopen(dist,5);
dist = imclose(dist,strel('sphere',1));
dist = imopen(dist,strel('sphere',1));
dist = imfill(dist,"holes");
dist = bwareaopen(dist,200);
cleaned(dist==1) = 5;
%% clean up proximal
prox = volTA==6;
prox = bwareaopen(prox,5);
prox = imclose(prox,strel('sphere',1));
prox = imopen(prox,strel('sphere',1));
prox = imfill(prox,"holes");
prox = bwareaopen(prox,200);
cleaned(prox==1) =6;
%% clean up henle
henle = volTA==7;
henle = bwareaopen(henle,5);
henle = imclose(henle,strel('sphere',2));
henle = imopen(henle,strel('sphere',1));
henle = imfill(henle,"holes");
henle = bwareaopen(henle,3000);
cleaned(henle==1) = 7;
%% clean up collecting duct
coll = volTA==8;
coll = bwareaopen(coll,5);
coll = imclose(coll,strel('sphere',2));
coll = imopen(coll,strel('sphere',1));
coll = imfill(coll,"holes");
coll = bwareaopen(coll,5000);
cleaned(coll==1) = 8;

%% clean up developing corpuscle
glom = volTA==10;
glom = imclose(glom,strel('sphere',2));
glom = imopen(glom,strel('sphere',1));
glom = bwareaopen(glom,250);
glom = imdilate(glom,strel('disk',2));
glom = imopen(glom,strel('disk',1));
cleaned(glom==1) = 10;

%% clean up undifferentiated blastema cells 
blas = volTA==11;
blas = imclose(blas,strel('sphere',2));
blas = imopen(blas,strel('sphere',1));
blas = bwareaopen(blas,10000);
blas = imdilate(blas,strel('disk',2));
blas = imopen(blas,strel('disk',1));
cleaned(blas==1) = 11;

%% clean up developing nephron 
dev = volTA==12;
dev = bwareaopen(dev,5);
dev = imclose(dev,strel('sphere',1));
dev = imopen(dev,strel('sphere',1));
dev = imfill(dev,"holes");
dev = bwareaopen(dev,200);
cleaned(dev==1) = 12;

%% clean up arteries
art = volTA==13;
art = bwareaopen(art,20);
art = imclose(art,strel('sphere',2));
art = imfill(art,'holes');
art = imdilate(art,strel('sphere',1));
art = imopen(art,strel('disk',1));
art = bwareaopen(art,2000);
cleaned(art==1) = 13;

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
save('\\Lucie Dequiedt\Kidney Project\Volumes\cleaned\E17_K1.mat','cleaned');