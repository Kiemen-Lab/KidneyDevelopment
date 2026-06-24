
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};

compositions = zeros(length(samples),4);
%% Load and pre-process
for j=1:length(samples)


    load([pth,samples{j},'.mat'])
    load([pth,'volzones2\',samples{j},'.mat'],'volzone');
    load(fullfile(pth, 'inside\', [samples{j} '.mat']),'volsurface');
    

    vol = volTA;
    clearvars -except vol volzone pth samples j compositions volsurface
    
    volzone = double(volzone);
    vol = double(vol);
    vol(vol==14) = 0;
    vol = vol.*double(volsurface);
    
     numzones = 3;
     for i=1:numzones
        zone = vol.*double(volzone==i);
        compositions(j,i) = sum(zone~=0,"all");
     end
    compositions(j,4) =sum(vol~=0,"all");
end

compo(1:length(samples)-4,1:3) = compositions(1:length(samples)-4,1:3)./compositions(1:length(samples)-4,4)*100;

compo(7,1:3) = sum(compositions(7:8,1:3),1)./sum(compositions(7:8,4),1)*100;
compo(8,1:3) = sum(compositions(9:10,1:3),1)./sum(compositions(9:10,4),1)*100;

