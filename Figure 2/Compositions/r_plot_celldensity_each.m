pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';

classNames = {'Urothelium','Med. Coll. Duct','Coll. Duct','Dist. Tub.','Henle Loop','Prox. Tub.', ...
              'Glom. Tuft','Bow. Cap.','Vein','Arteries','Arterioles','Stroma','Dev. Corpuscle','Dev. Nephron','Blastema'};

sample_labels = {'E17 (A)', 'E17 (B)', '15 weeks', '17 weeks', '19 wks (A)', '19 wks (B)', 'E80 (A)', 'E80 (B)'};

% Voxel volume in mm³: (4 µm)³ = 64 µm³ = 64 * 10^-9 mm³
voxel_volume_mm3 = 4 * 4 * 4 * 1e-9;  % mm³

%% Load saved data
load([pth, 'volume_and_cell_summary.mat'], 'volume_by_class', 'cells_by_class');
volume_by_class(:,14) = [];
cells_by_class(:,14) = [];
%% Compute cell density for each sample and tissue type (cells/mm³)
% volume_by_class: 8 samples x 16 tissue types (in voxels)
% cells_by_class: 8 samples x 16 tissue types (cell counts)

% Convert volume to mm³
volume_mm3 = volume_by_class * voxel_volume_mm3;

% Compute density (cells/mm³) for each sample and tissue
cell_density_matrix = zeros(size(volume_by_class));  % 8 samples x 15 tissue types
for i = 1:size(cell_density_matrix,1)
    for j = 1:size(cell_density_matrix,2)
        if volume_mm3(i, j) > 0
            cell_density_matrix(i, j) = cells_by_class(i, j) / volume_mm3(i, j);
        else
            cell_density_matrix(i, j) = NaN;  % No data
        end
    end
end

%% Heatmap 1: Cell Density by Sample and Tissue Type
figure('Position', [100 100 1200 700]);
imagesc(cell_density_matrix');
colormap(sky);
cb = colorbar;
ylabel(cb, 'Cell Density (cells/mm³)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:8, 'XTickLabel', sample_labels, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:15, 'YTickLabel', classNames);
xlabel('Sample', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Tissue Type', 'FontSize', 12, 'FontWeight', 'bold');
title('Cell Packing Density by Sample and Tissue Type', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);
axis tight;

% Add text values on heatmap
hold on;
for i = 1:size(cell_density_matrix,1)
    for j = 1:size(cell_density_matrix,2)
        if ~isnan(cell_density_matrix(i, j)) && cell_density_matrix(i, j) > 0
            text(i, j, sprintf('%.1e', cell_density_matrix(i, j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', '#a83232');
        end
    end
end
hold off;

%% Compute average cell density by species
% E17: samples 1-2
% Human: samples 3-6
% Macaque: samples 7-8

species_density = zeros(3, 15);  % 3 species x 15 tissue types
species_labels = {'Mouse', 'Human', 'Macaque'};

for j = 1:15
    % E17 (samples 1-2)
    e17_data = cell_density_matrix(1:2, j);
    e17_data = e17_data(~isnan(e17_data) & e17_data > 0);
    if ~isempty(e17_data)
        species_density(1, j) = mean(e17_data);
    else
        species_density(1, j) = NaN;
    end
    
    % Human (samples 3-6)
    human_data = cell_density_matrix(3:6, j);
    human_data = human_data(~isnan(human_data) & human_data > 0);
    if ~isempty(human_data)
        species_density(2, j) = mean(human_data);
    else
        species_density(2, j) = NaN;
    end
    
    % Macaque (samples 7-8)
    mac_data = cell_density_matrix(7:8, j);
    mac_data = mac_data(~isnan(mac_data) & mac_data > 0);
    if ~isempty(mac_data)
        species_density(3, j) = mean(mac_data);
    else
        species_density(3, j) = NaN;
    end
end

%% Heatmap 2: Average Cell Density by Species and Tissue Type
figure('Position', [100 100 800 700]);
imagesc(species_density');
colormap(sky);
cb = colorbar;
ylabel(cb, 'Cell Density (cells/mm³)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:3, 'XTickLabel', species_labels, 'XTickLabelRotation', 0);
set(gca, 'YTick', 1:15, 'YTickLabel', classNames);
xlabel('Species', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Tissue Type', 'FontSize', 12, 'FontWeight', 'bold');
title('Average Cell Packing Density by Species and Tissue Type', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);
axis tight;

% Add text values on heatmap
hold on;
for i = 1:3
    for j = 1:15
        if ~isnan(species_density(i, j)) && species_density(i, j) > 0
            text(i, j, sprintf('%.2e', species_density(i, j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', '#a83232', 'FontWeight', 'bold');
        end
    end
end
hold off;

%% Save results
