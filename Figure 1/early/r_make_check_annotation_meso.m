im = imread("\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 2\early\mouse\LD_MM11_chop_0503_he.tif");
J = imread("\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 2\early\mouse\LD_MM11_chop_0503_meso.tif");
mas = imread("\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 2\early\mouse\LD_MM11_chop_0503_urinary.tif");
mas = double(mas);
mas=mas==4;mas=imfill(mas,'holes');mas = imerode(mas,strel('disk',1));
J = double(J).*mas;%J(J==0) =5; 
ds=1;
cmap = [97 25 63;246 235 250;255 255 255];
  cmap2=cat(1,[0 0 0],cmap)/255;

im=im(1:ds:end,1:ds:end,:);J=J(1:ds:end,1:ds:end,:);
I=im2double(im);J=double(J);
J1=cmap2(J+1,1);J1=reshape(J1,size(J));
J2=cmap2(J+1,2);J2=reshape(J2,size(J));
J3=cmap2(J+1,3);J3=reshape(J3,size(J));
mask=cat(3,J1,J2,J3);
I2=(I*0.6)+(mask*0.4);
I2=uint8(I2*255);