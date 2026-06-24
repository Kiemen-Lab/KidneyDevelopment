%% Z-Projection Cell Density Heatmaps — Standalone Script
% Generates per-sample Z-projection heatmaps of glomerular cell density.
% Run this script directly; no function calls required.

clear; clc;

%% ── Configuration ────────────────────────────────────────────────────────────
pth     = '\\Lucie Dequiedt\Kidney Project\Volumes\';
outpth  = fullfile(pth, 'density_heatmaps');
if ~exist(outpth, 'dir'); mkdir(outpth); end

celldensityConversion = 4*4*4*1e-9;
distanceConversion    = 4/1000;
distanceBins          = 0:5:350;
smoothSigma           = 0;

volumeThresholds = struct( ...
    'E17_K1',        100000,  'E17_K2',        100000, ...
    'Mac_a',        1500000,  'Mac_b',        1500000, ...
    'Mac_c',        1500000,  'Mac_d',        1500000, ...
    'Hum_K1',       1500000,  'Hum_K2',       2000000, ...
    'Hum_K3_bottom',2500000,  'Hum_K3_top',   2500000);

% {sampleKey, displayLabel, species}
% Paired macaque kidneys are summed — listed with '+' separator
sampleList = { ...
    'E17_K1',        'E17 K1',     'mouse';
    'E17_K2',        'E17 K2',     'mouse';
    'Mac_a',   'Mac a',    'macaque';
    'Mac_b', 'Mac b', 'macaque';
    'Mac_c',   'Mac c',    'macaque';
    'Mac_d', 'Mac d', 'macaque';
    'Hum_K1',        'Hum K1',     'human';
    'Hum_K2',        'Hum K2',     'human';
    'Hum_K3_bottom', 'Hum K3 bot', 'human';
    'Hum_K3_top',    'Hum K3 top', 'human'};

speciesColors = struct( ...
    'mouse',   [217 189 213]/255, ...
    'macaque', [205  92  92]/255, ...
    'human',   [130  73  73]/255);

load(fullfile('\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 3\diameters_nuclei.mat'), 'D');
corrFactor = 4 / (4 + mean(D));

%% ── Sample loop ──────────────────────────────────────────────────────────────
nSamples = size(sampleList, 1);

for s = 1:nSamples

    key     = sampleList{s, 1};
    label   = sampleList{s, 2};
    species = sampleList{s, 3};
    fprintf('Processing %s...\n', label);

    %% Load data ───────────────────────────────────────────────────────────────
    if contains(key, '+')
        parts = strsplit(key, '+');
        sA = parts{1};  sB = parts{2};

        load(fullfile(pth, [sA '.mat']), 'volTA');
        mask = imread(fullfile(pth, 'distK_update\masks', [sA '.tif']));
        load(fullfile(pth, 'distK_update', [sA '.mat']), 'good');
        volTA  = double(volTA ~= 14) .* double(mask);
        distK  = bwdist(good);
        distK(volTA == 0) = 0;
        volTA  = imclose(volTA, strel('disk', 3));
        volTA  = imerode(volTA, strel('disk', 4));
        volTA  = bwareaopen(volTA, 100000);
        distK(volTA == 0) = 0;

        load(fullfile(pth, [sA '_cells.mat']), 'vc');
        vc = double(vc);
        load(fullfile(pth, [sB '_cells.mat']), 'vc2');
        vc = vc + double(vc2);

        thresh = volumeThresholds.(sA);

    else
        load(fullfile(pth, [key '.mat']), 'volTA');
        mask = imread(fullfile(pth, 'distK_update\masks', [key '.tif']));
        load(fullfile(pth, 'distK_update', [key '.mat']), 'good');
        volTA  = double(volTA ~= 14) .* double(mask);
        distK  = bwdist(good);
        distK(volTA == 0) = 0;

        if contains(key, {'K', 'M'})
            volTA = imclose(volTA, strel('disk', 3));
            volTA = imerode(volTA, strel('disk', 4));
            volTA = bwareaopen(volTA, 100000);
        end
        distK(volTA == 0) = 0;

        load(fullfile(pth, [key '_cells.mat']), 'vc');
        vc = double(vc);

        thresh = volumeThresholds.(key);
    end

    %% Build Z-projection density map ─────────────────────────────────────────
    [ny, nx, ~] = size(distK);
    countMap    = zeros(ny, nx);
    volumeMap   = zeros(ny, nx);

    for j = 1:length(distanceBins)-1
        minD = distanceBins(j);
        maxD = distanceBins(j+1);

        shellMask = (distK > minD) & (distK <= maxD);
        if sum(shellMask(:)) < thresh; continue; end
        cellMaps = vc.*shellMask;
        cellMaps = bwlabeln(cellMaps);
        countMap  = countMap  + max(cellMaps(:)) * corrFactor;
        volumeMap = volumeMap + sum(double(shellMask), 3);
    end

    densityMap             = NaN(ny, nx);
    validPix               = volumeMap > 0;
    densityMap(validPix)   = (countMap(validPix) ./ volumeMap(validPix)) / celldensityConversion;

    volTA_proj             = max(double(volTA > 0), [], 3);
    densityMap(volTA_proj == 0) = NaN;

    %% Smooth ──────────────────────────────────────────────────────────────────
    if smoothSigma > 0
        densMapSmooth = imgaussfilt(densityMap, smoothSigma, ...
                        'FilterDomain', 'spatial', 'Padding', 'replicate');
        densMapSmooth(isnan(densityMap)) = NaN;
    else
        densMapSmooth = densityMap;
    end

    %% Plot ────────────────────────────────────────────────────────────────────
    [ny, nx] = size(densMapSmooth);
    xmm = (1:nx) * distanceConversion;
    ymm = (1:ny) * distanceConversion;

    rgb    = speciesColors.(species);
    %cmap   = [linspace(1, rgb(1), 256)', ...
            %  linspace(1, rgb(2), 256)', ...
             % linspace(1, rgb(3), 256)'];

    fig = figure('Position', [100 100 800 650], 'Color', 'w', 'Visible', 'off');
    imagesc(xmm, ymm, densMapSmooth, 'AlphaData', ~isnan(densMapSmooth));
    axis image;
    colormap('turbo');
    set(gca, 'Color', [0.92 0.92 0.92], 'FontName', 'Arial', 'FontSize', 12);

    validVals = densMapSmooth(~isnan(densMapSmooth) & densMapSmooth > 0);
    if ~isempty(validVals)
        clim([prctile(validVals, 2), prctile(validVals, 98)]);
    end

    cb = colorbar;
    cb.Label.String  = 'Cell density [cells/mm³]';
    cb.Label.FontSize = 12;
    cb.FontSize       = 11;

    hold on;
    outline = boundarymask(double(~isnan(densityMap) & volTA_proj > 0));
    contour(xmm, ymm, outline, [0.5 0.5], 'k', 'LineWidth', 1);
    hold off;

    xlabel('X [mm]', 'FontSize', 13);
    ylabel('Y [mm]', 'FontSize', 13);
    title(sprintf('%s — Z-projection cell density', label), ...
          'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

    outFile = fullfile(outpth, [strrep(key, '+', '_') '_density_heatmap.png']);
    exportgraphics(fig, outFile, 'Resolution', 200);
    close(fig);
    fprintf('  Saved: %s\n', outFile);

    clearvars distK vc volTA volTA_proj countMap volumeMap densityMap densMapSmooth mask good;
end

%% ── Summary panel ────────────────────────────────────────────────────────────
fprintf('\nBuilding summary panel...\n');

cols = 4;
rows = ceil(nSamples / cols);
fig  = figure('Position', [100 100 cols*380 rows*340], 'Color', 'w', 'Visible', 'off');

for s = 1:nSamples

    key     = sampleList{s, 1};
    label   = sampleList{s, 2};
    species = sampleList{s, 3};

    if contains(key, '+')
        parts = strsplit(key, '+');
        sA = parts{1};  sB = parts{2};

        load(fullfile(pth, [sA '.mat']), 'volTA');
        mask  = imread(fullfile(pth, 'distK_update\masks', [sA '.tif']));
        load(fullfile(pth, 'distK_update', [sA '.mat']), 'good');
        volTA = double(volTA ~= 14) .* double(mask);
        distK = bwdist(good);
        distK(volTA == 0) = 0;
        volTA = imclose(volTA, strel('disk', 3));
        volTA = imerode(volTA, strel('disk', 4));
        volTA = bwareaopen(volTA, 100000);
        distK(volTA == 0) = 0;

        load(fullfile(pth, [sA '_cells.mat']), 'vc');
        vc = double(vc);
        load(fullfile(pth, [sB '_cells.mat']), 'vc2');
        vc = vc + double(vc2);
        thresh = volumeThresholds.(sA);

    else
        load(fullfile(pth, [key '.mat']), 'volTA');
        mask  = imread(fullfile(pth, 'distK_update\masks', [key '.tif']));
        load(fullfile(pth, 'distK_update', [key '.mat']), 'good');
        volTA = double(volTA ~= 14) .* double(mask);
        distK = bwdist(good);
        distK(volTA == 0) = 0;

        if contains(key, {'K', 'M'})
            volTA = imclose(volTA, strel('disk', 3));
            volTA = imerode(volTA, strel('disk', 4));
            volTA = bwareaopen(volTA, 100000);
        end
        distK(volTA == 0) = 0;

        load(fullfile(pth, [key '_cells.mat']), 'vc');
        vc     = double(vc);
        thresh = volumeThresholds.(key);
    end

    [ny, nx, ~] = size(distK);
    countMap    = zeros(ny, nx);
    volumeMap   = zeros(ny, nx);

    for j = 1:length(distanceBins)-1
        minD = distanceBins(j);
        maxD = distanceBins(j+1);
        shellMask = (distK > minD) & (distK <= maxD);
        if sum(shellMask(:)) < thresh; continue; end
        countMap  = countMap  + sum(double(vc > 0) .* double(shellMask), 3) * corrFactor;
        volumeMap = volumeMap + sum(double(shellMask), 3);
    end

    densityMap           = NaN(ny, nx);
    validPix             = volumeMap > 0;
    densityMap(validPix) = (countMap(validPix) ./ volumeMap(validPix)) / celldensityConversion;
    volTA_proj           = max(double(volTA > 0), [], 3);
    densityMap(volTA_proj == 0) = NaN;

    if smoothSigma > 0
        densMapSmooth = imgaussfilt(densityMap, smoothSigma, ...
                        'FilterDomain', 'spatial', 'Padding', 'replicate');
        densMapSmooth(isnan(densityMap)) = NaN;
    else
        densMapSmooth = densityMap;
    end

    [ny, nx] = size(densMapSmooth);
    xmm = (1:nx) * distanceConversion;
    ymm = (1:ny) * distanceConversion;

    rgb  = speciesColors.(species);
    cmap = [linspace(1, rgb(1), 256)', ...
            linspace(1, rgb(2), 256)', ...
            linspace(1, rgb(3), 256)'];

    ax = subplot(rows, cols, s);
    imagesc(ax, xmm, ymm, densMapSmooth, 'AlphaData', ~isnan(densMapSmooth));
    axis(ax, 'image');
    colormap(ax, cmap);
    set(ax, 'Color', [0.92 0.92 0.92], 'FontName', 'Arial', 'FontSize', 9);

    validVals = densMapSmooth(~isnan(densMapSmooth) & densMapSmooth > 0);
    if ~isempty(validVals)
        clim(ax, [prctile(validVals, 2), prctile(validVals, 98)]);
    end

    cb = colorbar(ax);
    cb.FontSize = 7;

    hold(ax, 'on');
    outline = boundarymask(double(~isnan(densityMap) & volTA_proj > 0));
    contour(ax, xmm, ymm, outline, [0.5 0.5], 'k', 'LineWidth', 0.8);
    hold(ax, 'off');

    title(ax, label, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    xlabel(ax, 'X [mm]', 'FontSize', 8);
    ylabel(ax, 'Y [mm]', 'FontSize', 8);

    clearvars distK vc volTA volTA_proj countMap volumeMap densityMap densMapSmooth mask good;
end

sgtitle('Z-Projection Cell Density — All Samples', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

outFile = fullfile(outpth, 'summary_panel_all_samples.png');
exportgraphics(fig, outFile, 'Resolution', 200);
close(fig);

fprintf('Done. All outputs saved to:\n  %s\n', outpth);