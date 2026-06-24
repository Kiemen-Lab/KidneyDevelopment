%% Paths and sample list
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

samples = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};

% Output directory for results
outpth = fullfile(pth, 'tubule_analysis');
if ~exist(outpth, 'dir')
    mkdir(outpth);
end

% Output directory for distK projections
distKpth = fullfile(pth, 'distK_tubules');
if ~exist(distKpth, 'dir')
    mkdir(distKpth);
end

%% Define distance bins (in voxels or your unit of measurement)
distanceBins = 0:5:350; % Adjust these values based on your kidney size

%% Loop over samples
for i = 1:length(samples)

    sampleName = samples{i};
    fprintf('Processing %s...\n', sampleName);
    
    % Create sample-specific directory for distK projections
    sampleDistKpth = fullfile(distKpth, sampleName);
    if ~exist(sampleDistKpth, 'dir')
        mkdir(sampleDistKpth);
    end
    
    %% Load sample-specific data
    % Load cleaned tubules volume (for tubule type identification)
    load(fullfile(pth, 'tubule_volumes_cleaned', [sampleName '.mat']));
    tubules = double(volTA);  % This contains cleaned tubule labels: 5,6,7,8
    
    % Load original volTA (for distance calculation and total volume)
    load(fullfile(pth, [sampleName '.mat']));
    volTA_original = double(volTA);  % Original volume for distance and total volume
       
    load(fullfile(pth, 'clean outside\', [sampleName '.mat']), 'outershell');
    load(fullfile(pth, 'inside\', [sampleName '.mat']), 'volsurface');
    clearvars good volTA
    
    %% Distance to kidney outer edge
    % Use original volTA for distance calculation (exclude label 14)
    volTA_bin = double(volTA_original ~= 14) .* double(volsurface);
    distK = bwdist(outershell);
    distK(volTA_bin == 0) = 0;
    
    if contains(sampleName,{'K';'M'})
           volTA_bin = imclose(volTA_bin,strel('disk',3));
           volTA_bin = imerode(volTA_bin,strel('disk',4));
           volTA_bin = bwareaopen(volTA_bin,100000);
    end
    distK(volTA_bin == 0) = 0;

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
    
    %% Compute tubule volume fractions at different distances
    nBins = length(distanceBins) - 1;
    volumes = zeros(nBins, 1);
    
    % Arrays to hold counts and fractions for each tubule type
    % Order: distal(5), proximal(6), loop(7), collecting duct(8)
    tubuleCounts = zeros(nBins, 4);
    tubuleFractions = zeros(nBins, 4);
    
    for j = 1:nBins
        % Define distance range
        minDist = distanceBins(j);
        maxDist = distanceBins(j+1);
        
        % Create mask for this distance shell using original volTA
        shellMask = (distK > minDist) & (distK <= maxDist);
        shellMask = shellMask & (volTA_bin > 0);  % ensure inside tissue
        
        %% Save z-projection for this distance bin
        shellMask_zprojMax = max(double(shellMask), [], 3);
        distK_masked = distK .* double(shellMask);
        distK_masked_zprojMax = max(distK_masked, [], 3);
        distK_masked_zprojMean = mean(distK_masked, 3);
        
        % Only save if there's actual data in this bin
        if sum(shellMask(:)) > 0
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
        end
        
        % Calculate total volume of tissue in this distance range (from original volTA)
        volumeVoxels = sum(shellMask(:));
        volumes(j) = volumeVoxels;
        
        % Compute tubule counts using CLEANED tubules volume
        if volumeVoxels > 0
            % Logical masks for each tubule type inside the shell
            % Use cleaned tubules for identification
            isDist = (tubules == 5) & shellMask;  % distal
            isProx = (tubules == 6) & shellMask;  % proximal
            isLoop = (tubules == 7) & shellMask;  % loop of Henle
            isColl = (tubules == 8) & shellMask;  % collecting duct
            
            cntDist = sum(isDist(:));
            cntProx = sum(isProx(:));
            cntLoop = sum(isLoop(:));
            cntColl = sum(isColl(:));
            
            tubuleCounts(j, :) = [cntDist, cntProx, cntLoop, cntColl];
            % Fraction = tubule voxels / total tissue voxels in shell
            tubuleFractions(j, :) = tubuleCounts(j, :) / volumeVoxels;
        else
            tubuleCounts(j, :) = 0;
            tubuleFractions(j, :) = 0;
        end
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
    
    %% Store results (only requested outputs)
    results.(sampleName).distanceBins = distanceBins;
    results.(sampleName).binCenters = distanceBins(1:end-1) + diff(distanceBins)/2;
    results.(sampleName).volumes = volumes;
    results.(sampleName).tubuleCounts = tubuleCounts;         % nBins x 4
    results.(sampleName).tubuleFractions = tubuleFractions;   % nBins x 4
    results.(sampleName).tubuleLabels = {'Distal(5)','Proximal(6)','Loop(7)','Collecting(8)'};
    
    %% Print summary statistics
    fprintf('  Z-projections saved to: %s\n', sampleDistKpth);
    
    % Print tubule summary (total voxels across all shells)
    totalTubCounts = sum(tubuleCounts, 1);
    totalTissueVox = sum(volumes);
    fprintf('  Total tissue volume: %d voxels\n', totalTissueVox);
    fprintf('  Tubule voxels (total): distal=%d, prox=%d, loop=%d, coll=%d\n', ...
        totalTubCounts(1), totalTubCounts(2), totalTubCounts(3), totalTubCounts(4));
    if totalTissueVox > 0
        fprintf('  Tubule fractions (global): distal=%.3f, prox=%.3f, loop=%.3f, coll=%.3f\n', ...
            totalTubCounts(1)/totalTissueVox, totalTubCounts(2)/totalTissueVox, ...
            totalTubCounts(3)/totalTissueVox, totalTubCounts(4)/totalTissueVox);
    end
    
    %% Plot for this sample
    figure('Name', sampleName, 'Position', [100 100 1400 600]);
    
    subplot(1,2,1);
    bar(results.(sampleName).binCenters, volumes);
    xlabel('Distance from surface (voxels)');
    ylabel('Volume (voxels)');
    title([sampleName ': Tissue volume per distance bin']);
    grid on;
    
    % Stacked bar of tubule volume fractions
    subplot(1,2,2);
    bf = results.(sampleName).tubuleFractions;
    h = bar(results.(sampleName).binCenters, bf, 'stacked');
    xlabel('Distance from surface (voxels)');
    ylabel('Volume fraction');
    title([sampleName ': Tubule volume fractions per distance bin']);
    legend(results.(sampleName).tubuleLabels, 'Location', 'bestoutside');
    grid on;
    
    % Save figure
    saveas(gcf, fullfile(outpth, [sampleName '_tubulefractions.png']));
    close(gcf);
    
    fprintf('  Completed %s\n\n', sampleName);
end

%% Save all results
save(fullfile(outpth, 'tubule_analysis_results.mat'), 'results');

%% Create summary plots for each tubule type across all samples
tubuleNames = {'Distal', 'Proximal', 'Loop of Henle', 'Collecting Duct'};

for tIdx = 1:4
    figure('Name', [tubuleNames{tIdx} ' Fraction Across Samples'], 'Position', [100 100 1200 800]);
    hold on;
    colors = lines(length(samples));
    for i = 1:length(samples)
        sn = samples{i};
        plot(results.(sn).binCenters, results.(sn).tubuleFractions(:, tIdx), ...
            '-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', colors(i,:), ...
            'DisplayName', sn);
    end
    xlabel('Distance from surface (voxels)');
    ylabel([tubuleNames{tIdx} ' volume fraction']);
    title([tubuleNames{tIdx} ' Tubule Volume Fraction - All Samples']);
    legend('Location', 'best');
    grid on;
    saveas(gcf, fullfile(outpth, ['all_samples_' lower(strrep(tubuleNames{tIdx}, ' ', '_')) '_fraction.png']));
    close(gcf);
end

%% Create grouped comparison for tubule fractions
for tIdx = 1:4
    figure('Name', [tubuleNames{tIdx} ' - Grouped Comparison'], 'Position', [100 100 1400 500]);
    
    % E17 samples
    subplot(1,3,1);
    hold on;
    for i = 1:length(samples)
        if contains(samples{i}, 'E17')
            plot(results.(samples{i}).binCenters, results.(samples{i}).tubuleFractions(:, tIdx), ...
                '-o', 'LineWidth', 1.5, 'DisplayName', samples{i});
        end
    end
    xlabel('Distance from surface (voxels)');
    ylabel([tubuleNames{tIdx} ' volume fraction']);
    title(['E17 Samples - ' tubuleNames{tIdx}]);
    legend('Location', 'best');
    grid on;
    
    % Human samples
    subplot(1,3,2);
    hold on;
    for i = 1:length(samples)
        if contains(samples{i}, 'Hum')
            plot(results.(samples{i}).binCenters, results.(samples{i}).tubuleFractions(:, tIdx), ...
                '-o', 'LineWidth', 1.5, 'DisplayName', samples{i});
        end
    end
    xlabel('Distance from surface (voxels)');
    ylabel([tubuleNames{tIdx} ' volume fraction']);
    title(['Human Samples - ' tubuleNames{tIdx}]);
    legend('Location', 'best');
    grid on;
    
    % Mac samples
    subplot(1,3,3);
    hold on;
    for i = 1:length(samples)
        if contains(samples{i}, 'Mac')
            plot(results.(samples{i}).binCenters, results.(samples{i}).tubuleFractions(:, tIdx), ...
                '-o', 'LineWidth', 1.5, 'DisplayName', samples{i});
        end
    end
    xlabel('Distance from surface (voxels)');
    ylabel([tubuleNames{tIdx} ' volume fraction']);
    title(['Mac Samples - ' tubuleNames{tIdx}]);
    legend('Location', 'best');
    grid on;
    
    saveas(gcf, fullfile(outpth, ['grouped_' lower(strrep(tubuleNames{tIdx}, ' ', '_')) '_comparison.png']));
    close(gcf);
end

%% Create combined tubule fraction plot for all samples (4 subplots)
figure('Name', 'All Tubule Fractions - All Samples', 'Position', [100 100 1400 1000]);

for tIdx = 1:4
    subplot(2,2,tIdx);
    hold on;
    colors = lines(length(samples));
    for i = 1:length(samples)
        sn = samples{i};
        plot(results.(sn).binCenters, results.(sn).tubuleFractions(:, tIdx), ...
            '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', colors(i,:), ...
            'DisplayName', sn);
    end
    xlabel('Distance from surface (voxels)');
    ylabel('Volume fraction');
    title(tubuleNames{tIdx});
    if tIdx == 1
        legend('Location', 'bestoutside');
    end
    grid on;
end

sgtitle('Tubule Volume Fractions - All Samples');
saveas(gcf, fullfile(outpth, 'all_tubules_all_samples.png'));

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results saved to: %s\n', outpth);
fprintf('Distance K projections saved to: %s\n', distKpth);
fprintf('Results structure saved to: %s\n', fullfile(outpth, 'tubule_analysis_results.mat'));