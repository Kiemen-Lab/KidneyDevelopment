pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};

%% Load and pre-process
for j=1:length(samples)
    % if exist([pth,'volzones\',samples{j},'.mat'])
    %     disp(['skip ',samples{j}]);
    %     continue
    % end
    S = load([pth,samples{j},'.mat']);
    vol = uint8(S.volTA);
    clear S

    volsurface = load(fullfile(pth, 'inside\', [samples{j} '.mat']),'volsurface');
    volsurface = volsurface.volsurface;
    % Build distance mask (logical)
    %dist = vol~=1 & vol~=9 & vol~=14;
    % dist = vol~=9 &vol~=14;
    % dist = imopen(dist, strel('sphere',3));
    % dist = imclose(dist, strel('sphere',5));
   
    
    % Apply mask in-place
    vol(~volsurface) = uint8(14);
    clear dist
    
    %% Build spherical kernel (single precision)
    rad0 = 150;
    rad  = round(rad0/4);
    
    dim = linspace(-1,1,rad*2+1);
    [X,Y,Z] = ndgrid(dim,dim,dim);
    sphere = double(sqrt(X.^2 + Y.^2 + Z.^2) <= 1);
    clear X Y Z dim
    
    volcircle = sum(sphere(:),'native');  % stays single
    
    %% Allocate output (single)
    volconv = zeros([size(vol), 3], 'single');
    
    %% Helper function logic inline (memory efficient)
    if contains(samples{j},'E17')
        labels = [15, 3, 12];  % medullary CD, glom, shell
    else
         labels = [8, 3, 12];  % CD, glom, shell
    end
    for i = 1:3
        mask = (vol == labels(i));   % logical
    
        mask = imclose(mask, strel('sphere',1));
        mask = imopen(mask,  strel('sphere',1));
        mask = bwareaopen(mask, 500);
    
        % Convolution (single)
        tmp = convn(single(mask), sphere, 'same');
        volconv(:,:,:,i) = tmp / volcircle * 100;
    
        clear mask tmp
    end
    
    %% Assign zones
    zz = all(volconv == 0, 4);
    
    [~, volzone] = max(volconv, [], 4);
    volzone = uint8(volzone);   % 1–3
    
    volzone(zz) = uint8(4);
    
    clear volconv zz

     save([pth,'volzones2\',samples{j},'.mat'],'volzone');
end
