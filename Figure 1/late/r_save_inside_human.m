pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};

for j=2%:length(samples)
   sampleName = samples{j};
    fprintf('Processing %s...\n', sampleName);
    S = load([pth,samples{j},'.mat']);
    vol = uint8(S.volTA);
    
    out = vol~=9 &vol~=14&vol~=2&vol~=13&vol~=16;
    i=1;
    if contains(sampleName,'Hum_K3_bottom')
        mask = imread(fullfile(pth, 'distK_update\masks', [samples{j} '.tif']));
        vol = double(vol).*mask;
        out = vol~=14&vol~=0;
        i=1;
    end

    if contains(sampleName,'Hum_K1')
        mask = imread(fullfile(pth, 'distK_update\masks', [samples{j} '.tif']));
        vol = double(vol).*mask;
        out = vol~=14&vol~=0;
        i=2;
    end

     if contains(sampleName,'Hum_K3_top')
        out = vol~=14&vol~=0;
        i=1;
     end

     if contains(sampleName,'Mac')
         mask = imread(fullfile(pth, 'distK_update\masks', [samples{j} '.tif']));
        vol = double(vol).*mask;
         out = vol~=0&vol~=14&vol~=1;
         i=1;
    end
   
    volsurface=uint8(zeros(size(out)));
    st = [100000 5000];
    
    for z=1:size(out,3)
       tmp=out(:,:,z);
       tmp = tmp~=0;
       tmp = bwareaopen(tmp,1000);
       tmp=imclose(tmp,strel('disk',10));
       tmp=imfill(tmp,'holes');
       tmp = bwareaopen(tmp,st(i));
       volsurface(:,:,z)=tmp;
       disp([z size(out,3)]);
      
    end 
    
     save(fullfile(pth,'inside\', [sampleName '.mat']),'volsurface','-v7.3')
end
