%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { 'Hum_K2'...
   'Hum_K3_bottom' 'Hum_K3_top'};

% Output directory for distK projections
distKpth = fullfile(pth, 'clean outside');


for i = 1:length(samples)
    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);
    if exist(fullfile(distKpth, [sampleName '.mat']),'file')
        disp('already processed')
        continue
    end
    load([pth,'volzones2\',sampleName,'.mat']);
    load([pth,sampleName,'.mat']);
    load([pth,'inside\',sampleName,'.mat'])
   
   out = volTA==12|volTA==11|volTA==10;

outershell = bwperim(volsurface,26);

 surface=uint8(zeros(size(out)));

    for z=1:size(out,3)
       tmp=out(:,:,z);
       % tmp = imopen(tmp,strel('disk',5));
        tmp = bwareaopen(tmp,100);
       tmp = imclose(tmp,strel('disk',30));
    %  tmp = imfill(tmp,'holes');
        tmp = bwareaopen(tmp,1000);

       surface(:,:,z)=tmp;
       disp([z size(out,3)]);
      
    end 

 surface = bwdist(~surface);
 surface = surface>15&surface<20;
 outershell = surface==1|outershell==1;


 outershell(:,:,1) = outershell(:,:,2);
 save(fullfile(distKpth, [sampleName '.mat']),'outershell','-v7.3')
end
