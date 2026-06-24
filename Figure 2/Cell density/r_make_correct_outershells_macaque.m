%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};


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
    out = out==11|out==12|out==10;

    volsurface=uint8(zeros(size(out)));
    figure;
    for z=1:size(out,3)
       tmp=out(:,:,z);
       tmp = tmp~=0;
       tmp=imclose(tmp,strel('disk',30));
       tmp=imfill(tmp,'holes');
       tmp = bwareafilt(tmp,1);
       volsurface(:,:,z)=tmp;
       disp([z size(out,3)]);
      
    end 

     a=1;
    if  contains(sampleName,'_a')
        for z=1:79
        volsurface(:,:,z) = volsurface(:,:,80);
        end
        for z=87:91
        volsurface(:,:,z) = volsurface(:,:,86);

        end
         volsurface(:,:,103) = volsurface(:,:,104);
         volsurface(:,:,126) = volsurface(:,:,125);
    end
     if  contains(sampleName,'_b')
        for z=1:64
        volsurface(:,:,z) = volsurface(:,:,65);
        end
        volsurface(:,:,69) = volsurface(:,:,68);
     end
      if  contains(sampleName,'_c')
        for z=1:68
        volsurface(:,:,z) = volTA(:,:,z)~=14&volTA(:,:,z)~=9&volTA(:,:,z)~=1;
        volsurface(:,:,z) =  imclose(volsurface(:,:,z),strel('disk',30));
        volsurface(:,:,z)=imfill(volsurface(:,:,z),'holes');
        end
      end

      if  contains(sampleName,'_d')
        for z=1:15
        volsurface(:,:,z) = volsurface(:,:,16);
        end
        volsurface(:,:,21) = volsurface(:,:,22);
         volsurface(:,:,50) = volsurface(:,:,51);
        volsurface(:,:,63) = volsurface(:,:,64);
        volsurface(:,:,80) = volsurface(:,:,81);
        volsurface(:,:,201) = volsurface(:,:,203);
        volsurface(:,:,202) = volsurface(:,:,203);
        volsurface(:,:,205) = volsurface(:,:,209);
        volsurface(:,:,208) = volsurface(:,:,206);
        volsurface(:,:,214) = volsurface(:,:,213);
     end

   
    t = bwperim(volsurface,26);
    t = bwareaopen(t,500);
    outershell = t;

    outershell(:,:,1) = outershell(:,:,2);
    save(fullfile(distKpth, [sampleName '.mat']),'outershell','-v7.3')
end
