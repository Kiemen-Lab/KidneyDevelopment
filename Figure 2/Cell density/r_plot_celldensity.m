%% Load and Setup
load('\\Lucie Dequiedt\Kidney Project\Volumes\celldensity\cell_density_analysis.mat')

% Configuration
config.distanceConversion = 4/1000; % Convert to mm
config.volumeThresholds = struct(...
    'E17_K1', 100000, 'E17_K2', 100000, ...
    'Mac_a', 1500000, 'Mac_b', 1500000, ...
    'Mac_c', 1500000, 'Mac_d', 1500000, ...
    'Hum_K1', 1500000, 'Hum_K2', 2000000, ...
    'Hum_K3_bottom', 2500000, 'Hum_K3_top', 2500000);
config.celldensityConversion = 4*4*4*10^(-9);

% Color scheme by species
colors = struct(...
    'mouse', [217, 189, 213]/255, ...
    'macaque', [205, 92, 92]/255, ...
    'human', [130, 73, 73]/255);

% Age information - define as separate variables
age_mouse = 'E17';
age_macaque = '80 days';
age_human = {'15 weeks', '17 weeks', '19 weeks', '19 weeks'}; % K1, K2, K3_bottom, K3_top

%% Process All Samples
plotData = struct();

% Mouse samples
[dist1, dens1] = processData(results, 'E17_K1', config.volumeThresholds.E17_K1, config.distanceConversion,config.celldensityConversion);
plotData.mouse{1} = struct('name', 'E17 K1', 'dist', dist1, 'dens', dens1);

[dist2, dens2] = processData(results, 'E17_K2', config.volumeThresholds.E17_K2, config.distanceConversion,config.celldensityConversion);
plotData.mouse{2} = struct('name', 'E17 K2', 'dist', dist2, 'dens', dens2);

% Macaque samples (combined pairs)
% Mac_a + Mac_b
distances = results.Mac_a.binCenters(:)' * config.distanceConversion;
volumes_ab = results.Mac_a.volumes + results.Mac_b.volumes;
counts_ab = results.Mac_a.correctedCounts + results.Mac_b.correctedCounts;
dens_ab = counts_ab ./ (volumes_ab*config.celldensityConversion);
validIdx = volumes_ab >= config.volumeThresholds.Mac_a;
plotData.macaque{1} = struct('name', 'Mac a+b', 'dist', distances(validIdx), 'dens', dens_ab(validIdx));

% Mac_c + Mac_d
volumes_cd = results.Mac_c.volumes + results.Mac_d.volumes;
counts_cd = results.Mac_c.correctedCounts + results.Mac_d.correctedCounts;
dens_cd = counts_cd ./ (volumes_cd*config.celldensityConversion);
validIdx = volumes_cd >= config.volumeThresholds.Mac_c;
plotData.macaque{2} = struct('name', 'Mac c+d', 'dist', distances(validIdx), 'dens', dens_cd(validIdx));

% Human samples
humanSamples = {'Hum_K1', 'Hum_K2', 'Hum_K3_bottom', 'Hum_K3_top'};
humanNames = {'Hum K1', 'Hum K2', 'Hum K3 bottom', 'Hum K3 top'};
plotData.human = cell(1, length(humanSamples));
for i = 1:length(humanSamples)
    [dist_h, dens_h] = processData(results, humanSamples{i}, ...
        config.volumeThresholds.(humanSamples{i}), config.distanceConversion,config.celldensityConversion);
    plotData.human{i} = struct('name', humanNames{i}, 'dist', dist_h, 'dens', dens_h);
end

figure('Position', [100, 100, 900, 600], 'Color', 'w');
hold on; box on; 

% Plot mouse data
for i = 1:length(plotData.mouse)
    if i == 1
        plot(plotData.mouse{i}.dist, plotData.mouse{i}.dens, ...
            '-o', 'Color', colors.mouse, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors.mouse, ...
            'DisplayName', 'Mouse');
    else
        plot(plotData.mouse{i}.dist, plotData.mouse{i}.dens, ...
            '-o', 'Color', colors.mouse, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors.mouse, ...
            'HandleVisibility', 'off');
    end
end

% Plot macaque data
for i = 1:length(plotData.macaque)
    if i == 1
        plot(plotData.macaque{i}.dist, plotData.macaque{i}.dens, ...
            '-o', 'Color', colors.macaque, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors.macaque, ...
            'DisplayName', 'Macaque');
    else
        plot(plotData.macaque{i}.dist, plotData.macaque{i}.dens, ...
            '-o', 'Color', colors.macaque, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors.macaque, ...
            'HandleVisibility', 'off');
    end
end

% Plot human data
for i = 1:length(plotData.human)
    if i == 1
        plot(plotData.human{i}.dist, plotData.human{i}.dens, ...
            '-o', 'Color', colors.human, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors.human, ...
            'DisplayName', 'Human');
    else
        plot(plotData.human{i}.dist, plotData.human{i}.dens, ...
            '-o', 'Color', colors.human, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors.human, ...
            'HandleVisibility', 'off');
    end
end

% Formatting
xlabel('Distance to Nephrogenic Zone [mm]', 'FontSize', 15);
ylabel('Cell Density [cells/mm^3]', 'FontSize', 15);
legend('Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 15, 'LineWidth', 1.5);
set(gca, 'FontName', 'Arial')

% Add age labels next to curves
% Mouse age label (place near first mouse curve)
if ~isempty(plotData.mouse)
    maxDist_mouse = plotData.mouse{1}.dist(end);
    maxDens_mouse = 1.05*plotData.mouse{1}.dens(end);
    text(maxDist_mouse, maxDens_mouse, ['  ' age_mouse], ...
        'Color', colors.mouse, 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle');
end

% Macaque age label (place near first macaque curve)
if ~isempty(plotData.macaque)
    maxDist_mac = plotData.macaque{1}.dist(end-6);
    maxDens_mac = 1.3*plotData.macaque{1}.dens(end);
    text(maxDist_mac, maxDens_mac, ['  ' age_macaque], ...
        'Color', colors.macaque, 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle');
end

% Human age labels (one for each sample)
for i = 1:length(plotData.human)
    if ~isempty(plotData.human{i}.dist)
        maxDist_hum = plotData.human{i}.dist(end);
        maxDens_hum = plotData.human{i}.dens(end);
        text(maxDist_hum, maxDens_hum, ['  ' age_human{i}], ...
            'Color', colors.human, 'FontSize', 12, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'middle');
    end
end

%% Save plots (optional)
% saveas(gcf, 'cell_density_comparison.png');
% saveas(gcf, 'cell_density_comparison.fig');

%% Local Functions (must be at the end)
function [dist_filtered, dens_filtered] = processData(results, sampleName, threshold, distanceConversion,celldensityConversion)
    % Extract distances
    distances = results.(sampleName).binCenters(:)' * distanceConversion;
    
    % Get cell density (or calculate from counts/volumes)
    if isfield(results.(sampleName), 'cellDensity')
        cellDensity = results.(sampleName).cellDensity(:)';
    else
        cellDensity = results.(sampleName).correctedCounts(:)' ./ results.(sampleName).volumes(:)';
    end
    
    % Get volumes and apply threshold
    volumes = results.(sampleName).volumes(:)';
    validIdx = volumes >= threshold;
    
    dist_filtered = distances(validIdx);
    dens_filtered = cellDensity(validIdx)/celldensityConversion;
end