%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { ...
    'Hum_K1'};

% Output directory for distK projections
distKpth = fullfile(pth, 'clean outside');

%% Output directory for distK projections
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
 

   out = double(volTA);
   out = out~=14;
   

    volsurface=uint8(zeros(size(out)));

    for z=1:size(out,3)
       tmp=out(:,:,z);
       tmp = imopen(tmp,strel('disk',5));
       tmp = bwareaopen(tmp,10000);
       tmp = imclose(tmp,strel('disk',5));
       tmp = bwareaopen(tmp,10000);
       tmp = imfill(tmp,'holes');
       volsurface(:,:,z)=tmp;
       disp([z size(out,3)]);
      
    end 

   
    t = bwperim(volsurface,26);
    t = bwareaopen(t,500);
    outershell = t;

    outershell(:,:,1) = outershell(:,:,2);
    save(fullfile(distKpth, [sampleName '.mat']),'outershell','-v7.3')
end
