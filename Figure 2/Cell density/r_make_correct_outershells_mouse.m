%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { ...
    'E17_K1', 'E17_K2'};

% Output directory for glomeruli distributions
outpth = fullfile(pth, 'glomeruli');
if ~exist(outpth, 'dir')
    mkdir(outpth);
end

% Output directory for distK projections
distKpth = fullfile(pth, 'distK_update\masks');
if ~exist(distKpth, 'dir')
    mkdir(distKpth);
end

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


    out = double(volTA).*double(volzone==3);
    out(out==14 | out==1 | out==9) =0;
    dist = bwdist(out~=0);

    volsurface=uint8(zeros(size(out)));
    figure;
    for z=1:size(out,3)
       tmp=out(:,:,z);
       tmp = tmp~=0;
       tmp=imclose(tmp,strel('disk',10));
       tmp=imfill(tmp,'holes');
       volsurface(:,:,z)=tmp;
       imshow(tmp);
       disp([z size(out,3)]);
    end 

   
    t = bwperim(volsurface,26);
    figure;
    for z=1:size(out,3)
          imshow(t(:,:,z));
    end
    t = bwareaopen(t,200);
    outershell = t;

        save(fullfile(distKpth, [sampleName '.mat']),'outershell','-v7.3')
end
