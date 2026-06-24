pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};

for j=1:length(samples)
   sampleName = samples{j};
    fprintf('Processing %s...\n', sampleName);
    S = load([pth,samples{j},'.mat']);
    vol = uint8(S.volTA);
    
    out = vol~=1 &vol~=9 &vol~=14&vol~=2&vol~=13&vol~=16;
    volsurface=uint8(zeros(size(out)));
    
    for z=1:size(out,3)
       tmp=out(:,:,z);
       tmp = tmp~=0;
       tmp=imclose(tmp,strel('disk',5));
       tmp=imfill(tmp,'holes');
       tmp = bwareaopen(tmp,1000);
       volsurface(:,:,z)=tmp;
       disp([z size(out,3)]);
      
    end 
    
     save(fullfile(pth,'inside\', [sampleName '.mat']),'volsurface','-v7.3')
end
