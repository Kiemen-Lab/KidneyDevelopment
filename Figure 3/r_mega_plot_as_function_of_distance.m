%% Paths
pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
datapth = fullfile(pth, 'glomeruli');
densitypth = fullfile(pth, 'density_surface');

%% Distance binning (in mm)
binWidth_mm = 0.05;
edges_mm = 0:binWidth_mm:1.5;
binCenters_mm = edges_mm(1:end-1) + binWidth_mm/2;

minGlomPerBin = 30;

%% Species definition
species(1).name = 'Mouse';
species(1).pattern = 'E17';

species(2).name = 'Human';
species(2).pattern = 'Hum';

species(3).name = 'Macaque';
species(3).pattern = 'Mac';

cols = [217, 189, 213; ...
    130, 73, 73; ...
    205, 92, 92]/255; 

%% Define human sample labels
humanLabels = struct(...
    'Hum_K1', '15 weeks', ...
    'Hum_K2', '17 weeks', ...
    'Hum_K3_top', '19 weeks (A)', ...
    'Hum_K3_bottom', '19 weeks (B)');

%% ========== PART 1: Load Glomerular Volume Data ==========
files = dir(fullfile(datapth, '*_mature_glomeruli_distributions.mat'));
samples_vol = struct;
sampleIdx = 0;

for f = 1:length(files)
    fname = files(f).name;
    load(fullfile(datapth, fname), 'distances', 'volumes');
    
    %% Identify species
    sp = [];
    for s = 1:length(species)
        if contains(fname, species(s).pattern)
            sp = s;
            break
        end
    end
    if isempty(sp); continue; end
    
    %% Handle Macaque samples
    if sp == 3
        if contains(fname, 'Mac_a', 'IgnoreCase', true)
            mac_a_distances = distances;
            mac_a_volumes = volumes;
            continue;
        elseif contains(fname, 'Mac_b', 'IgnoreCase', true)
            distances = [mac_a_distances; distances];
            volumes = [mac_a_volumes; volumes];
            sampleName = 'Mac 1';
        elseif contains(fname, 'Mac_c', 'IgnoreCase', true)
            mac_c_distances = distances;
            mac_c_volumes = volumes;
            continue;
        elseif contains(fname, 'Mac_d', 'IgnoreCase', true)
            distances = [mac_c_distances; distances];
            volumes = [mac_c_volumes; volumes];
            sampleName = 'Mac 2';
        else
            continue;
        end
    else
        sampleName = erase(fname,'_mature_glomeruli_distributions.mat');
    end
    
    dist_mm = (distances * 4) / 1000;
    volumes = volumes / 1e9;
    
    %% Bin volumes
    meanVol = nan(length(binCenters_mm),1);
    stdVol  = nan(length(binCenters_mm),1);
    
    for b = 1:length(binCenters_mm)
        inBin = dist_mm >= edges_mm(b) & dist_mm < edges_mm(b+1);
        v = volumes(inBin);
        if numel(v) >= minGlomPerBin
            meanVol(b) = mean(v);
            stdVol(b)  = std(v)/ sqrt(numel(v));
        end
    end
    
    sampleIdx = sampleIdx + 1;
    samples_vol(sampleIdx).name = sampleName;
    samples_vol(sampleIdx).species = sp;
    samples_vol(sampleIdx).meanVol = meanVol;
    samples_vol(sampleIdx).stdVol  = stdVol;
    samples_vol(sampleIdx).dist_mm = binCenters_mm;
end

%% ========== PART 2: Load Cell Density Data ==========
files = dir(fullfile(datapth, '*_mature_glomeruli_distributions.mat'));
samples_cell = struct;
sampleIdx = 0;

for f = 1:length(files)
    fname = files(f).name;
    load(fullfile(datapth, fname), 'distances', 'volumes', 'cellcount');
    
    %% Identify species
    sp = [];
    for s = 1:length(species)
        if contains(fname, species(s).pattern)
            sp = s;
            break
        end
    end
    if isempty(sp); continue; end
    
    %% Handle Macaque samples
    if sp == 3
        if contains(fname, 'Mac_a', 'IgnoreCase', true)
            mac_a_distances = distances;
            mac_a_volumes = volumes;
            mac_a_cellcount = cellcount;
            continue;
        elseif contains(fname, 'Mac_b', 'IgnoreCase', true)
            distances = [mac_a_distances; distances];
            volumes = [mac_a_volumes; volumes];
            cellcount = [mac_a_cellcount; cellcount];
            sampleName = 'Mac 1';
        elseif contains(fname, 'Mac_c', 'IgnoreCase', true)
            mac_c_distances = distances;
            mac_c_volumes = volumes;
            mac_c_cellcount = cellcount;
            continue;
        elseif contains(fname, 'Mac_d', 'IgnoreCase', true)
            distances = [mac_c_distances; distances];
            volumes = [mac_c_volumes; volumes];
            cellcount = [mac_c_cellcount; cellcount];
            sampleName = 'Mac 2';
        else
            continue;
        end
    else
        sampleName = erase(fname,'_mature_glomeruli_distributions.mat');
    end
    
    dist_mm = (distances * 4) / 1000;
    volumes_mm3 = volumes / 1e9;
    celldens = cellcount ./ volumes_mm3;
    
    %% Bin cell density
    meanDens = nan(length(binCenters_mm),1);
    semDens  = nan(length(binCenters_mm),1);
    
    for b = 1:length(binCenters_mm)
        inBin = dist_mm >= edges_mm(b) & dist_mm < edges_mm(b+1);
        v = celldens(inBin);
        if numel(v) >= minGlomPerBin
            meanDens(b) = mean(v);
            semDens(b)  = std(v) / sqrt(numel(v));
        end
    end
    
    sampleIdx = sampleIdx + 1;
    samples_cell(sampleIdx).name = sampleName;
    samples_cell(sampleIdx).species = sp;
    samples_cell(sampleIdx).meanDens = meanDens;
    samples_cell(sampleIdx).semDens  = semDens;
    samples_cell(sampleIdx).dist_mm = binCenters_mm;
end

%% ========== PART 3: Load Density Surface Data ==========
samples_density = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b', 'Mac_c', 'Mac_d'};

bincenters = {};
vols = {};
vessel_vf = {};
gloms = {};

for i = 1:length(samples_density)
    sampleName = samples_density{i};
    respth = fullfile(densitypth, sampleName);
    load([respth,'\density_analysis.mat']);
    
    bincenters{i} = results_table.Distance_Min + (results_table.Distance_Max - results_table.Distance_Min)/2;
    vols{i} = shell_volumes;
    vessel_vf{i} = vessel_volume_fractions;
    gloms{i} = glom_density;
end

%% Combine macaque samples
arta = vols{7}.*vessel_vf{7}+vols{8}.*vessel_vf{8};
gloma = vols{7}.*gloms{7}+vols{8}.*gloms{8};
vola = vols{7} +vols{8};
vessel_vf{7} = arta./vola;
gloms{7} = gloma./vola;
vols{7} = vola;
bincenters{7} = bincenters{7}; % Keep first one

artb = vols{9}.*vessel_vf{9}+vols{10}.*vessel_vf{10};
glomb = vols{9}.*gloms{9}+vols{10}.*gloms{10};
volb = vols{9} +vols{10};
vessel_vf{8} = artb./volb;
gloms{8} = glomb./volb;
vols{8} = volb;
bincenters{8} = bincenters{9}; % Keep first one

gloms(9:10) = [];
vessel_vf(9:10) = [];
vols(9:10) = [];
bincenters(9:10) = [];

config.volumeThresholds = struct(...
    'E17_K1', 0.1, 'E17_K2', 0.1, ...
    'Mac_a', 0.4, 'Mac_b', 0.4, ...
    'Hum_K1', 0.2, 'Hum_K2', 0.2, ...
    'Hum_K3_bottom', 0.4, 'Hum_K3_top', 0.4);

% Keep original order but we'll plot in different order
samples_plot = { ...
    'E17_K1', 'E17_K2', ...
    'Hum_K1', 'Hum_K2', ...
    'Hum_K3_bottom', 'Hum_K3_top', ...
    'Mac_a', 'Mac_b'};

% Create plot order indices: Macaque (7,8), Human (3,4,5,6), Mouse (1,2)
plot_order_indices = [7, 8, 3, 4, 5, 6, 1, 2];

%% ========== CREATE COMBINED FIGURE (NO NORMALIZATION) ==========
figure('Position', [100, 100, 800, 1200]);

plotOrder = [3 2 1];   % Macaque → Human → Mouse

%% SUBPLOT 1: Glomerular Volume
subplot(4, 1, 1);
hold on;

legendHandles = gobjects(length(species),1);

for s = plotOrder
    idx = find([samples_vol.species] == s);
    
    for k = idx
        % Use absolute distance (no normalization)
        dist_valid = samples_vol(k).dist_mm(~isnan(samples_vol(k).meanVol));
        vol_valid = samples_vol(k).meanVol(~isnan(samples_vol(k).meanVol));
        std_valid = samples_vol(k).stdVol(~isnan(samples_vol(k).meanVol));
        
        if ~isempty(dist_valid)
            errorbar(dist_valid, vol_valid, std_valid, ...
                '-o', 'LineWidth', 3, 'Color', cols(s,:), ...
                'HandleVisibility', 'off', 'MarkerSize', 3, 'MarkerFaceColor', cols(s,:));
            
            % % Add label for human samples
            % if s == 2 && isfield(humanLabels, samples_vol(k).name)
            %     text(dist_valid(end) + 0.02, vol_valid(end), ...
            %         humanLabels.(samples_vol(k).name), ...
            %         'FontSize', 9, 'Color', cols(s,:), ...
            %         'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            % end
        end
    end
    
    legendHandles(s) = plot(nan, nan, '-o', ...
        'Color', cols(s,:), 'LineWidth', 3, ...
        'DisplayName', species(s).name, 'MarkerSize', 3, 'MarkerFaceColor', cols(s,:));
end

ylabel('Glom. Volume [mm³]', 'FontSize', 10);
legend(legendHandles, 'Location', 'best', 'FontSize', 10);
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'XTickLabel', []);
box on;
hold off;

%% SUBPLOT 2: Cell Density
subplot(4, 1, 2);
hold on;
%add = 0.02;
for s = plotOrder
    idx = find([samples_cell.species] == s);
    
    for k = idx
        % Use absolute distance (no normalization)
        dist_valid = samples_cell(k).dist_mm(~isnan(samples_cell(k).meanDens));
        dens_valid = samples_cell(k).meanDens(~isnan(samples_cell(k).meanDens));
        sem_valid = samples_cell(k).semDens(~isnan(samples_cell(k).meanDens));
        
        if ~isempty(dist_valid)
            errorbar(dist_valid, dens_valid, sem_valid, ...
                '-o', 'LineWidth', 3, 'Color', cols(s,:), ...
                'HandleVisibility', 'off', 'MarkerSize', 3, 'MarkerFaceColor', cols(s,:));
            
            % Add label for human samples
            if s == 2 && isfield(humanLabels, samples_cell(k).name)
                if contains(samples_cell(k).name,'K2')
                    add = 0.02;
                    yy = 0.9;
                else 
                    add = 0.02;
                    yy = 1;
                end
                text(dist_valid(end) + add, yy*dens_valid(end), ...
                    humanLabels.(samples_cell(k).name), ...
                    'FontSize', 9, 'Color', cols(s,:), ...
                    'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            end
        end
    end
end

ylabel('Glom. Cell Dens. \newline [cells/mm³]', 'FontSize', 10);
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'XTickLabel', []);
box on;
hold off;

%% SUBPLOT 3: Glomeruli Density
subplot(4, 1, 3);
hold on;

% Plot in order: Macaque → Human → Mouse
for i = plot_order_indices
    sampleName = samples_plot{i};
    
    if startsWith(sampleName, 'E17')
        color = cols(1,:);
    elseif startsWith(sampleName, 'Mac')
        color = cols(3,:);
    elseif startsWith(sampleName, 'Hum')
        color = cols(2,:);
    end
    
    threshold = config.volumeThresholds.(sampleName);
    valid_idx = vols{i} >= threshold;
    
    % Use absolute distance in mm (no normalization)
    distance_mm = bincenters{i}*4/1000;
    
    plot(distance_mm(valid_idx), gloms{i}(valid_idx), ...
        'o-', 'Color', color, 'LineWidth', 3, 'MarkerSize', 3, ...
        'MarkerFaceColor', color, 'HandleVisibility', 'off');
    
    
end

ylabel('Glom. Density \newline [glom/mm³]', 'FontSize', 10);
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'XTickLabel', []);
box on;
hold off;

%% SUBPLOT 4: Vessel Volume Fraction
subplot(4, 1, 4);
hold on;

% Plot in order: Macaque → Human → Mouse
for i = plot_order_indices
    sampleName = samples_plot{i};
    
    if startsWith(sampleName, 'E17')
        color = cols(1,:);
    elseif startsWith(sampleName, 'Mac')
        color = cols(3,:);
    elseif startsWith(sampleName, 'Hum')
        color = cols(2,:);
    end
    
    threshold = config.volumeThresholds.(sampleName);
    valid_idx = vols{i} >= threshold;
    
    % Use absolute distance in mm (no normalization)
    distance_mm = bincenters{i}*(4/1000);
    
    plot(distance_mm(valid_idx), vessel_vf{i}(valid_idx)*100, ...
        'o-', 'Color', color, 'LineWidth', 3, 'MarkerSize', 3, ...
        'MarkerFaceColor', color, 'HandleVisibility', 'off');
    
end

xlabel('Distance from Nephrogenic Zone [mm]', 'FontSize', 13);
ylabel('Vessel V. F. [%]', 'FontSize', 10);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
box on;
hold off;