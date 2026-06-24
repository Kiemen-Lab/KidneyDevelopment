%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
pthglom = '\\Lucie Dequiedt\Kidney Project\Data for paper\Codes\Figure 4\Glom Volumes\volgloms\';
samples = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};

% Output directory for glomeruli distributions
outpth = fullfile(pth, 'density_surface');
if ~exist(outpth, 'dir')
    mkdir(outpth);
end

% Distance parameters
distance_step = 5;  % voxels between distance bins
max_distance = 350;  % maximum distance to analyze
distance_bins = 5:distance_step:max_distance;

% Voxel size for volume calculations (adjust to your actual voxel size)
voxel_size = [4, 4, 4];  % [x, y, z] in micrometers
voxel_volume = prod(voxel_size);  % volume of one voxel in µm³

%% Loop over samples
for i = 1:2%length(samples)
    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);
              
    % Create sample-specific output directory
    sample_outpth = fullfile(outpth, sampleName);
    if ~exist(sample_outpth, 'dir')
        mkdir(sample_outpth);
    end
    
    % Load sample-specific data
    load(fullfile(pthglom, [sampleName '_labelled_glom.mat']));
    load(fullfile(pth, [sampleName '.mat']));          % contains volTA
    load(fullfile(pth, 'clean outside\', [sampleName '.mat']), 'outershell');
    load(fullfile(pth, 'inside\', [sampleName '.mat']), 'volsurface');
    %% Distance to kidney outer edge
    v = double(volTA~=14).*double(volsurface);
    dist = bwdist(outershell);
    dist(v==0) = 0;
    
    if contains(sampleName,{'K';'M'})
           v = imclose(v,strel('disk',3));
           v = imerode(v,strel('disk',4));
           v = bwareaopen(v,100000);
          % distK = bwdist(~volTA);
    end
    dist(v==0) =0;
    clearvars v 
    vess = volTA==13|volTA==16;
    vess = bwareaopen(vess,25);

    sz = size(volTA);
    
    % Prepare vessel z-projection for visualization
    vess_zproj = max(vess, [], 3);
    
    % Initialize results arrays
    num_bins = length(distance_bins) - 1;
    glom_counts = zeros(num_bins, 1);
    glom_density = zeros(num_bins, 1);  % glomeruli per mm³
    glom_volume_fraction = zeros(num_bins, 1);  % fraction of space occupied by glomeruli
    vessel_volume_fractions = zeros(num_bins, 1);
    shell_volumes = zeros(num_bins, 1);  % volume of each shell
    
    % Loop over distance bins
    for j = 1:num_bins
        dist_min = distance_bins(j);
        dist_max = distance_bins(j+1);
        
        fprintf('  Distance %d to %d voxels...\n', dist_min, dist_max);
        
        % Create mask for current distance shell
        shell_mask = (dist >= dist_min) & (dist < dist_max);
        
        % Calculate shell volume
        total_voxels_in_shell = sum(shell_mask(:));
        shell_volume_mm3 = (total_voxels_in_shell * voxel_volume) / 1e9;  % Convert µm³ to mm³
        shell_volumes(j) = shell_volume_mm3;
        
        % Count glomeruli in this shell (MODIFIED TO ALLOW DOUBLE COUNTING)
        % Instead of counting unique glomerulus IDs, count all glomerulus voxels
        glom_in_shell = volglom .* shell_mask;
        glom_voxels_in_shell = sum(glom_in_shell(:) > 0);
        
        % For backwards compatibility, also count unique glomeruli
        glom_ids = unique(glom_in_shell);
        glom_ids(glom_ids == 0) = [];  % Remove background
        glom_counts(j) = length(glom_ids);
        
        % Calculate glomeruli density based on VOXELS (allows double counting)
        % This measures glomerular tissue density rather than object count
        if shell_volume_mm3 > 0
            % Density based on voxel volume rather than object count
            glom_voxel_volume_mm3 = (glom_voxels_in_shell * voxel_volume) / 1e9;
            glom_density(j) = glom_voxel_volume_mm3 / shell_volume_mm3 * 1000;  % Normalized per mm³
        else
            glom_density(j) = 0;
        end
        
        % Calculate glomeruli volume fraction
        if total_voxels_in_shell > 0
            glom_volume_fraction(j) = glom_voxels_in_shell / total_voxels_in_shell;
        else
            glom_volume_fraction(j) = 0;
        end
        
        % Calculate vessel volume fraction in this shell
        vess_in_shell = vess & shell_mask;
        vessel_voxels_in_shell = sum(vess_in_shell(:));
        
        if total_voxels_in_shell > 0
            vessel_volume_fractions(j) = vessel_voxels_in_shell / total_voxels_in_shell;
        else
            vessel_volume_fractions(j) = 0;
        end
        
        % Create visualization
        % Z-projections
        shell_zproj = max(shell_mask, [], 3);
        glom_shell_zproj = max(glom_in_shell > 0, [], 3);
        outershell_zproj = max(outershell, [], 3);
        
        % Create RGB composite image
        rgb_img = zeros(sz(1), sz(2), 3);
        rgb_img(:,:,1) = vess_zproj;              % Vessels in red
        rgb_img(:,:,2) = double(glom_shell_zproj); % Glomeruli in shell in green
        rgb_img(:,:,3) = double(shell_zproj) * 0.5; % Current shell in blue
        
        % Add outershell boundary in yellow (red + green)
        outershell_boundary = imdilate(outershell_zproj, strel('disk', 2)) & ~outershell_zproj;
        rgb_img(:,:,1) = max(rgb_img(:,:,1), double(outershell_boundary));
        rgb_img(:,:,2) = max(rgb_img(:,:,2), double(outershell_boundary));
        
        % Normalize and save
        if max(rgb_img(:)) > 0
            rgb_img = rgb_img / max(rgb_img(:));
        end
        
        % Create figure
        fig = figure('Visible', 'off');
        imshow(rgb_img);
        title(sprintf('%s: Distance %d-%d voxels from surface\nGlom Density: %.2f/mm³, Count: %d, VesselFrac: %.3f', ...
            strrep(sampleName, '_', '\_'), dist_min, dist_max, ...
            glom_density(j), glom_counts(j), vessel_volume_fractions(j)));
        
        % Save figure
        saveas(fig, fullfile(sample_outpth, sprintf('dist_%03d_%03d.png', dist_min, dist_max)));
        close(fig);
    end
    
    % Save results
    results_table = table(distance_bins(1:end-1)', distance_bins(2:end)', ...
        glom_counts, glom_density, glom_volume_fraction, vessel_volume_fractions, shell_volumes, ...
        'VariableNames', {'Distance_Min', 'Distance_Max', 'Glomeruli_Count', ...
        'Glomeruli_Density_per_mm3', 'Glomeruli_Volume_Fraction', 'Vessel_Volume_Fraction', 'Shell_Volume_mm3'});
    
    writetable(results_table, fullfile(sample_outpth, 'density_analysis.csv'));
    
    % Save mat file with results
    save(fullfile(sample_outpth, 'density_analysis.mat'), ...
        'distance_bins', 'glom_counts', 'glom_density', 'glom_volume_fraction', ...
        'vessel_volume_fractions', 'shell_volumes', 'results_table');
    
    % Create summary plot
    fig = figure('Visible', 'off', 'Position', [100, 100, 800, 900]);
    
    subplot(3,1,1);
    bar(distance_bins(1:end-1) + distance_step/2, glom_density);
    xlabel('Distance from kidney surface (voxels)');
    ylabel('Glomeruli density (per mm³)');
    title(sprintf('%s: Glomeruli density distribution', strrep(sampleName, '_', '\_')));
    grid on;
    
    subplot(3,1,2);
    bar(distance_bins(1:end-1) + distance_step/2, glom_volume_fraction);
    xlabel('Distance from kidney surface (voxels)');
    ylabel('Glomeruli volume fraction');
    title('Glomeruli volume fraction');
    grid on;
    
    subplot(3,1,3);
    bar(distance_bins(1:end-1) + distance_step/2, vessel_volume_fractions);
    xlabel('Distance from kidney surface (voxels)');
    ylabel('Vessel volume fraction');
    title('Vessel density');
    grid on;
    
    saveas(fig, fullfile(sample_outpth, 'summary_plots.png'));
    close(fig);
    
    fprintf('  Completed %s\n', sampleName);
end

fprintf('All samples processed!\n');