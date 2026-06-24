%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};

% Output directory for glomeruli distributions
outpth = fullfile(pth, 'glomeruli');
if ~exist(outpth, 'dir')
    mkdir(outpth);
end

%% Load constants (shared across samples)
load("diameters_nuclei.mat", 'D');
 
ses = [2 2 2 1 1 1 2 2 2 2];
%% Loop over samples
for i = 1:length(samples)

    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);
    % if ~exist(fullfile(outpth, [sampleName '_mature_glomeruli_distributions.mat']))
    %% Load sample-specific data
    load(fullfile(pth, [sampleName '.mat']));          % contains volTA
    load(fullfile(pth, [sampleName '_cells.mat']));    % contains vc
    load(fullfile(pth, 'outershells', [sampleName '.mat']), 'outershell');

    %% Distance to kidney outer edge
    distK = bwdist(outershell);

    %% Extract and clean glomeruli volume
    % Glomeruli are labeled as 3 or 10 in volTA
    % volglom = (volTA == 3) | (volTA == 10) | (volTA==4);
    volglom = (volTA == 3) | (volTA==10);

   
    se2 = strel('sphere', 2);
    volglom = imopen(volglom, se2);
    volglom = bwareaopen(volglom, 100);
    % Morphological cleanup
    se1 = strel('sphere', 2);
    volglom = imclose(volglom, se1);


    % Remove small objects (noise)
    volglom = bwareaopen(volglom, 500);
    volglom = imfill(volglom,'holes');

    % add bowmans capsule 
   volglom = bwlabeln(volglom);  % each glom is labeled 
  

    save(fullfile('\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 4\Glom Volumes\volgloms\', [sampleName '_labelled_glom.mat']), ...
        'volglom');
    %% Compute glomerulus metrics
    % [distances, volumes, cellcount, indices] = ...
    %     get_dist_centroid(volglom, distK, vc, D, [3 4]);
    % 
    % %% Save results for this sample
    % save(fullfile(outpth, [sampleName '_mature_glomeruli_distributions.mat']), ...
    %      'distances', 'volumes', 'cellcount', 'indices');
    % end
end
