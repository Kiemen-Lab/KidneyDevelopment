pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

% Output directory for results
outpth = fullfile(pth, 'tubule_volumes_cleaned');
if ~exist(outpth, 'dir')
    mkdir(outpth);
end

samples = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};

for i=1:length(samples)
    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);

    % Load sample-specific data
    load(fullfile(pth, [sampleName '.mat']));          % contains volTA
    dist = volTA==5;

        d = bwareaopen(dist,10);
        d = imclose(d,strel('sphere',3));
        d = imopen(d,strel('sphere',1));
        d = bwareaopen(d,50);
        
        prox = volTA==6;
        p = bwareaopen(prox,10);
        p = imclose(p,strel('sphere',3));
        p = imopen(p,strel('sphere',1));
        p = bwareaopen(p,50);
        
        
        hen = volTA==7;
        h = bwareaopen(hen,50);
        h = imclose(h,strel('sphere',3));
        h = imopen(h,strel('sphere',1));
        h = bwareaopen(h,100);
        
        coll = volTA==8;
        c = bwareaopen(coll,50);
        c = imclose(c,strel('sphere',3));
        c = imopen(c,strel('sphere',1));
        c = bwareaopen(c,100);

    volTA = zeros(size(volTA));
    volTA(d==1) = 5;
    volTA(p==1) = 6;
    volTA(h==1) = 7;
    volTA(c==1) = 8;

    save(fullfile(outpth, [sampleName '.mat']), ...
        'volTA', '-v7.3');

end