%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};

% Output directory for glomeruli distributions
outpth = fullfile(pth, 'celldensity');
if ~exist(outpth, 'dir')
    mkdir(outpth);
end

% Output directory for distK projections
distKpth = fullfile(pth, 'distK_new');
if ~exist(distKpth, 'dir')
    mkdir(distKpth);
end

%% Load constants (shared across samples)
load("diameters_nuclei.mat", 'D');

%% Define distance bins (in voxels or your unit of measurement)
distanceBins = 0:5:350; % Adjust these values based on your kidney size
% e.g., 0-10, 10-20, 20-30, ... voxels from surface

%% Loop over samples
for i = [1:length(samples)]

    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);
    
    % Create sample-specific directory for distK projections
    sampleDistKpth = fullfile(distKpth, sampleName);
    if ~exist(sampleDistKpth, 'dir')
        mkdir(sampleDistKpth);
    end
    
    %% Load sample-specific data
    load(fullfile(pth, [sampleName '.mat']));          % contains volTA
        
    load(fullfile(pth, [sampleName '_cells.mat']));    % contains vc
    vc = double(vc);
    load(fullfile(pth, 'inside\', [sampleName '.mat']));
    load(fullfile(pth, 'clean outside\', [sampleName '.mat']));
    
    clearvars good
    %% Distance to kidney outer edge
    volTA = double(volTA~=14).*double(volsurface);
    distK = bwdist(outershell);
    distK(volTA==0) = 0;
    
    if contains(sampleName,{'K';'M'})
           volTA = imclose(volTA,strel('disk',3));
           volTA = imerode(volTA,strel('disk',4));
           volTA = bwareaopen(volTA,100000);
          % distK = bwdist(~volTA);
    end
    distK(volTA==0) =0;

    %% Save full distK z-projection
    distK_zprojMax = max(distK, [], 3);
    distK_zprojMean = mean(distK, 3);
    
    figure('Visible', 'off');
    subplot(1,2,1);
    imagesc(distK_zprojMax);
    colorbar;
    title([sampleName ': distK - Max Z-projection']);
    axis equal tight;
    
    subplot(1,2,2);
    imagesc(distK_zprojMean);
    colorbar;
    title([sampleName ': distK - Mean Z-projection']);
    axis equal tight;
    
    saveas(gcf, fullfile(sampleDistKpth, [sampleName '_distK_full_zproj.png']));
    close(gcf);
    
    %% Get distance values at cell locations
    % Since vc is binary, we can directly index where vc==1
    cellDistances = distK(vc > 0);
    
    %% Compute cell density at different distances
    nBins = length(distanceBins) - 1;
    cellDensity = zeros(nBins, 1);
    rawCounts = zeros(nBins, 1);
    correctedCounts = zeros(nBins, 1);
    volumes = zeros(nBins, 1);
    
    for j = 1:nBins
        % Define distance range
        minDist = distanceBins(j);
        maxDist = distanceBins(j+1);
        
        % Create mask for this distance shell
        shellMask = (distK > minDist) & (distK <= maxDist);
        
        % Calculate volume of tissue in this distance range
        volumeVoxels = sum(shellMask(:));
        volumes(j) = volumeVoxels;
        
        % Skip computation if shellMask is empty
        if volumeVoxels == 0
            fprintf('  Bin %d (%d-%d voxels): Empty shell, skipping...\n', j, minDist, maxDist);
            cellDensity(j) = 0;
            rawCounts(j) = 0;
            correctedCounts(j) = 0;
            continue; % Skip to next iteration
        end
        
        %% Save z-projection for this distance bin (only if non-empty)
        shellMask_zprojMax = max(double(shellMask), [], 3);
        distK_masked = distK .* double(shellMask);
        distK_masked_zprojMax = max(distK_masked, [], 3);
        distK_masked_zprojMean = mean(distK_masked, 3);
        
        figure('Visible', 'off', 'Position', [100 100 1400 400]);
        
        % Shell mask projection
        subplot(1,3,1);
        imagesc(shellMask_zprojMax);
        colorbar;
        title(sprintf('Shell Mask: %d-%d voxels', minDist, maxDist));
        axis equal tight;
        colormap(gca, 'gray');
        
        % Max projection of distK in this shell
        subplot(1,3,2);
        imagesc(distK_masked_zprojMax);
        colorbar;
        title(sprintf('distK Max: %d-%d voxels', minDist, maxDist));
        axis equal tight;
        colormap(gca, 'jet');
        caxis([minDist, maxDist]);
        
        % Mean projection of distK in this shell
        subplot(1,3,3);
        imagesc(distK_masked_zprojMean);
        colorbar;
        title(sprintf('distK Mean: %d-%d voxels', minDist, maxDist));
        axis equal tight;
        colormap(gca, 'jet');
        
        % Save with distance range in filename
        filename = sprintf('%s_distK_shell_%04d_%04d.png', sampleName, minDist, maxDist);
        saveas(gcf, fullfile(sampleDistKpth, filename));
        close(gcf);
        
        % Count cells in this distance range (only if shell is not empty)
        rawCount = vc.*shellMask;
        rawCount = sum(rawCount(:));
        rawCounts(j) = rawCount;
        
        % Apply correction factor
        correctedCount = rawCount * (4 / (4 + mean(D)));
        correctedCounts(j) = correctedCount;
        
        % Calculate density (cells per unit volume)
        cellDensity(j) = correctedCount / volumeVoxels;
    end
    
    %% Create composite figure showing all distance shells
    % Find non-empty bins
    nonEmptyBins = find(volumes > 0);
    if ~isempty(nonEmptyBins)
        nPanels = min(length(nonEmptyBins), 20); % Limit to 20 panels
        step = max(1, floor(length(nonEmptyBins) / nPanels));
        selectedBins = nonEmptyBins(1:step:end);
        
        nRows = ceil(sqrt(length(selectedBins)));
        nCols = ceil(length(selectedBins) / nRows);
        
        figure('Visible', 'off', 'Position', [100 100 300*nCols 300*nRows]);
        for k = 1:length(selectedBins)
            j = selectedBins(k);
            minDist = distanceBins(j);
            maxDist = distanceBins(j+1);
            
            shellMask = (distK >= minDist) & (distK < maxDist);
            shellMask_zproj = max(double(shellMask), [], 3);
            
            subplot(nRows, nCols, k);
            imagesc(shellMask_zproj);
            title(sprintf('%d-%d', minDist, maxDist));
            axis equal tight off;
            colormap(gca, 'gray');
        end
        sgtitle([sampleName ': Distance Shells Overview']);
        saveas(gcf, fullfile(sampleDistKpth, [sampleName '_distK_shells_overview.png']));
        close(gcf);
    end
    
    %% Store results
    results.(sampleName).distanceBins = distanceBins;
    results.(sampleName).binCenters = distanceBins(1:end-1) + diff(distanceBins)/2;
    results.(sampleName).rawCounts = rawCounts;
    results.(sampleName).correctedCounts = correctedCounts;
    results.(sampleName).volumes = volumes;
    results.(sampleName).cellDensity = cellDensity;
    
    %% Print summary statistics
    fprintf('  Total cells (raw): %d\n', sum(vc(:)));
    fprintf('  Total cells (corrected): %.1f\n', sum(correctedCounts));
    fprintf('  Mean density (non-empty bins): %.6f cells/voxel\n', mean(cellDensity(volumes > 0)));
    fprintf('  Number of non-empty bins: %d / %d\n', sum(volumes > 0), nBins);
    fprintf('  Z-projections saved to: %s\n', sampleDistKpth);
    
    %% Plot for this sample
    figure('Name', sampleName, 'Position', [100 100 1200 800]);
    
    subplot(2,2,1);
    bar(results.(sampleName).binCenters, rawCounts);
    xlabel('Distance from surface (voxels)');
    ylabel('Raw cell count');
    title([sampleName ': Raw counts']);
    grid on;
    
    subplot(2,2,2);
    bar(results.(sampleName).binCenters, correctedCounts);
    xlabel('Distance from surface (voxels)');
    ylabel('Corrected cell count');
    title([sampleName ': Corrected counts']);
    grid on;
    
    subplot(2,2,3);
    bar(results.(sampleName).binCenters, volumes);
    xlabel('Distance from surface (voxels)');
    ylabel('Volume (voxels)');
    title([sampleName ': Tissue volume']);
    grid on;
    
    subplot(2,2,4);
    % Only plot non-zero density values
    validIdx = cellDensity > 0;
    plot(results.(sampleName).binCenters(validIdx), cellDensity(validIdx), ...
        '-o', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Distance from surface (voxels)');
    ylabel('Cell density (cells/voxel)');
    title([sampleName ': Cell density profile']);
    grid on;
    
    % Save figure
    saveas(gcf, fullfile(outpth, [sampleName '_density_profile.png']));
    
    fprintf('  Completed %s\n\n', sampleName);
end

%% Save all results
save(fullfile(outpth, 'cell_density_analysis.mat'), 'results');

%% Create summary plot across all samples
figure('Name', 'All Samples Comparison', 'Position', [100 100 1200 800]);
hold on;
legendEntries = {};
colors = lines(length(samples));
for i = 1:length(samples)
    sampleName = samples{i};
    % Only plot non-zero density values
    validIdx = results.(sampleName).cellDensity > 0;
    plot(results.(sampleName).binCenters(validIdx), ...
        results.(sampleName).cellDensity(validIdx), ...
        '-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', colors(i,:), ...
        'DisplayName', sampleName);
end
xlabel('Distance from kidney surface (voxels)');
ylabel('Cell density (cells/voxel)');
title('Cell Density Profiles - All Samples');
legend('Location', 'best');
grid on;
saveas(gcf, fullfile(outpth, 'all_samples_density_comparison.png'));

%% Create grouped comparison by sample type
figure('Name', 'Grouped Comparison', 'Position', [100 100 1400 500]);

% E17 samples
subplot(1,3,1);
hold on;
for i = 1:length(samples)
    if contains(samples{i}, 'E17')
        validIdx = results.(samples{i}).cellDensity > 0;
        plot(results.(samples{i}).binCenters(validIdx), ...
            results.(samples{i}).cellDensity(validIdx), ...
            '-o', 'LineWidth', 1.5, 'DisplayName', samples{i});
    end
end
xlabel('Distance from surface (voxels)');
ylabel('Cell density (cells/voxel)');
title('E17 Samples');
legend('Location', 'best');
grid on;

% Human samples
subplot(1,3,2);
hold on;
for i = 1:length(samples)
    if contains(samples{i}, 'Hum')
        validIdx = results.(samples{i}).cellDensity > 0;
        plot(results.(samples{i}).binCenters(validIdx), ...
            results.(samples{i}).cellDensity(validIdx), ...
            '-o', 'LineWidth', 1.5, 'DisplayName', samples{i});
    end
end
xlabel('Distance from surface (voxels)');
ylabel('Cell density (cells/voxel)');
title('Human Samples');
legend('Location', 'best');
grid on;

% Mac samples
subplot(1,3,3);
hold on;
for i = 1:length(samples)
    if contains(samples{i}, 'Mac')
        validIdx = results.(samples{i}).cellDensity > 0;
        plot(results.(samples{i}).binCenters(validIdx), ...
            results.(samples{i}).cellDensity(validIdx), ...
            '-o', 'LineWidth', 1.5, 'DisplayName', samples{i});
    end
end
xlabel('Distance from surface (voxels)');
ylabel('Cell density (cells/voxel)');
title('Mac Samples');
legend('Location', 'best');
grid on;

saveas(gcf, fullfile(outpth, 'grouped_density_comparison.png'));

fprintf('\nAnalysis complete! Results saved to: %s\n', outpth);
fprintf('Distance K projections saved to: %s\n', distKpth);