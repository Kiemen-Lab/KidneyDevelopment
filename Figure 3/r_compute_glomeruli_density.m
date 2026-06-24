%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
pthglom = '\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 4\Glom Volumes\volgloms\';
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
density = {};
%% Loop over samples
for i = 1:length(samples)

    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);
    % Load sample-specific data
    load(fullfile(pthglom, [sampleName '_labelled_glom.mat']));
    load(fullfile(pth, [sampleName '.mat']));          % contains volTA
    load([pth,'inside\',sampleName,'.mat'])
   
   
   
   numglom = max(volglom(:));
   volTA = double(volTA).*double(volsurface);
   volTA(volTA==0) = 14;
   voltot = sum(volTA(:)~=14);
   im = zeros(size(volglom,1:2));
   density{i,1} = numglom;
   density{i,2} = voltot*4*4*4*10^-9;

end
 save(fullfile(outpth, ['glomeruli_density.mat']), ...
         'density');


%% Load density data
load(fullfile(outpth, 'glomeruli_density.mat'), 'density');

%% Define sample groupings
species_map = {
    'E17', [1, 2];                    % E17_K1, E17_K2
    'Hum', [3, 4, 5, 6];              % Hum_K1, K2, K3_bottom, K3_top
    'Mac', [7, 8, 9, 10]              % Mac_a, b, c, d
};

%% Compute density per species

for i=1:length(samples)-4
    species = samples{i};
    total_gloms = cell2mat(density(i,1));
    total_volume = cell2mat(density(i,2));
    glom_density = total_gloms/total_volume;
    fprintf('%s: %d gloms / %.3f mm³ = %.2f gloms/mm³\n', ...
        species, total_gloms, total_volume, glom_density);
end


 species = 'Mac_a';
total_gloms = sum(cell2mat(density(7:8,1)));
    total_volume = sum(cell2mat(density(8:9,2)));
    glom_density = total_gloms/total_volume;
    fprintf('%s: %d gloms / %.3f mm³ = %.2f gloms/mm³\n', ...
        species, total_gloms, total_volume, glom_density);


     species = 'Mac_b';
total_gloms = sum(cell2mat(density(9:10,1)));
    total_volume = sum(cell2mat(density(9:10,2)));
    glom_density = total_gloms/total_volume;
    fprintf('%s: %d gloms / %.3f mm³ = %.2f gloms/mm³\n', ...
        species, total_gloms, total_volume, glom_density);