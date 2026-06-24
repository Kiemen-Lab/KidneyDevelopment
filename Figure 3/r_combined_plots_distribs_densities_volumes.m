pth = '\\Lucie Dequiedt\Kidney Project\Volumes\glomeruli\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};

%% Load data
for i=1:length(samples)
    load([pth,samples{i},'_mature_glomeruli_distributions.mat']);
    volumes_all{i} = volumes / 1e9;  % Convert to mm3
    celldens = cellcount./(volumes/1e9);
    celldens_all{i} = celldens;
end

% Define color schemes for each species
mouse_cols = [245, 223, 223; 194, 169, 191]/255;
hum_cols = [171, 169, 169; 168, 131, 131; 130, 73, 73; 56, 34, 34]/255;
mac_cols = [139, 0, 0; 178, 34, 34; 205, 92, 92; 240, 128, 128]/255;

% Shared limits for volumes
xlim_vol = [0, 7e-4];
ylim_vol = [0 0.3];
edges_vol = linspace(xlim_vol(1), xlim_vol(2), 25);

% Shared limits for cell density
xlim_dens = [0, 2.5e6];
ylim_dens = [0 0.6];
edges_dens = linspace(xlim_dens(1), xlim_dens(2), 25);

%% Create figure with 2 columns and 3 rows
figure('Position', [100, 100, 1600, 900]);

% Column widths and positions
col1_left = 0.08;
col2_left = 0.54;
col_width = 0.38;
row_height = 0.23;
row3_bottom = 0.11;
row2_bottom = 0.405;
row1_bottom = 0.70;

%% LEFT COLUMN - VOLUMES

% Plot 1: Mouse volumes (top left)
subplot('Position', [col1_left, row1_bottom, col_width, row_height]);
hold on
histogram(volumes_all{1}, edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.6, 'EdgeColor',mouse_cols(1,:),'FaceColor',mouse_cols(1,:));
histogram(volumes_all{2}, edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.1, 'EdgeColor',mouse_cols(2,:),'FaceColor',mouse_cols(1,:));
set(gca, 'FontName','Arial', 'FontSize',15, 'Box','on', 'LineWidth',1, 'XTickLabel', []);
legend({'E17 (A)','E17 (B)'}, 'Location','best');
title('Mouse');
xlim(xlim_vol);
ylim(ylim_vol);

% Plot 2: Human volumes (middle left)
subplot('Position', [col1_left, row2_bottom, col_width, row_height]);
hold on
histogram(volumes_all{3}, edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(1,:),'FaceColor',hum_cols(1,:));
histogram(volumes_all{4}, edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(2,:),'FaceColor',hum_cols(1,:));
histogram(volumes_all{5}, edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(3,:),'FaceColor',hum_cols(1,:));
histogram(volumes_all{6}, edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(4,:),'FaceColor',hum_cols(1,:));
set(gca, 'FontName','Arial', 'FontSize',15, 'Box','on', 'LineWidth',1, 'XTickLabel', []);
ylabel('Probability Distribution');
legend({'15 wks','17 wks','19 wks (A)','19 wks (B)'}, 'Location','best');
title('Human');
xlim(xlim_vol);
ylim(ylim_vol);

% Plot 3: Macaque volumes (bottom left)
subplot('Position', [col1_left, row3_bottom, col_width, row_height]);
hold on
histogram([volumes_all{7}.',volumes_all{8}.'], edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',mac_cols(4,:),'FaceColor',mac_cols(4,:));
histogram([volumes_all{9}.',volumes_all{10}.'], edges_vol, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',mac_cols(2,:),'FaceColor',mac_cols(4,:));
set(gca, 'FontName','Arial', 'FontSize',15, 'Box','on', 'LineWidth',1);
xlabel('Glomerular tuft volume [mm^3]');
legend({'E80 (A)','E80 (B)'}, 'Location','best');
title('Macaque');
xlim(xlim_vol);
ylim(ylim_vol);

%% RIGHT COLUMN - CELL DENSITY

% Plot 4: Mouse cell density (top right)
subplot('Position', [col2_left, row1_bottom, col_width, row_height]);
hold on
histogram(celldens_all{1}, edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.6, 'EdgeColor',mouse_cols(1,:),'FaceColor',mouse_cols(1,:));
histogram(celldens_all{2}, edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.1, 'EdgeColor',mouse_cols(2,:),'FaceColor',mouse_cols(1,:));
set(gca, 'FontName','Arial', 'FontSize',15, 'Box','on', 'LineWidth',1, 'XTickLabel', []);
legend({'E17 (A)','E17 (B)'}, 'Location','best');
title('Mouse');
xlim(xlim_dens);
ylim(ylim_dens);

% Plot 5: Human cell density (middle right)
subplot('Position', [col2_left, row2_bottom, col_width, row_height]);
hold on
histogram(celldens_all{3}, edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(1,:),'FaceColor',hum_cols(1,:));
histogram(celldens_all{4}, edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(2,:),'FaceColor',hum_cols(1,:));
histogram(celldens_all{5}, edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(3,:),'FaceColor',hum_cols(1,:));
histogram(celldens_all{6}, edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',hum_cols(4,:),'FaceColor',hum_cols(1,:));
set(gca, 'FontName','Arial', 'FontSize',15, 'Box','on', 'LineWidth',1, 'XTickLabel', []);
legend({'15 wks','17 wks','19 wks (A)','19 wks (B)'}, 'Location','best');
title('Human');
xlim(xlim_dens);
ylim(ylim_dens);

% Plot 6: Macaque cell density (bottom right)
subplot('Position', [col2_left, row3_bottom, col_width, row_height]);
hold on
histogram([celldens_all{7}.',celldens_all{8}.'], edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',mac_cols(4,:),'FaceColor',mac_cols(4,:));
histogram([celldens_all{9}.',celldens_all{10}.'], edges_dens, ...
    'Normalization','probability', ...
    'FaceAlpha',0.2, 'EdgeColor',mac_cols(2,:),'FaceColor',mac_cols(4,:));
set(gca, 'FontName','Arial', 'FontSize',15, 'Box','on', 'LineWidth',1);
xlabel('Cell Density [cells/mm^3]');
legend({'E80 (A)','E80 (B)'}, 'Location','best');
title('Macaque');
xlim(xlim_dens);
ylim(ylim_dens);