%% Comprehensive comparison script for E11 and MA40 tables
load("\\Lucie Dequiedt\Kidney Project\Volumes\metanephros\E11_meta.mat")
E11 = information;
load("\\Lucie Dequiedt\Kidney Project\Volumes\metanephros\MA40_meta.mat")
MA40 = information;

%% Define color scheme
% Base colors for each compartment
cmap_base = [215, 191, 155;...   % 1  ureteric_tree
             135 76 57;...     % 2  renal_vesicle
             193 174 167;...     % 3  condensed_mesenchyme (cap)
             242, 228, 231]/255;   % 4  loose_stroma

% Create lighter shade for mouse (E11) and darker shade for macaque (MA40)
cmap_mouse = min(cmap_base + 0.15, 1);  % Lighter
cmap_macaque = max(cmap_base - 0.15, 0); % Darker

%% Voxel size conversion
voxel_size_um3 = 4 * 4 * 4;  % um^3
voxel_to_mm3 = voxel_size_um3 / (1000^3);  % Convert to mm^3

%% Extract data for E11
% K1
E11_K1_TotVol = E11{1, 'TotVol'}{1} * voxel_to_mm3;
E11_K1_TreeVol = E11{1, 'TreeVol'}{1} * voxel_to_mm3;
E11_K1_CapVol = E11{1, 'CapVol'}{1} * voxel_to_mm3;
E11_K1_StrVol = E11{1, 'StrVol'}{1} * voxel_to_mm3;
E11_K1_TipsVol = sum(E11{1, 'TipsVolumes'}{1}) * voxel_to_mm3;
E11_K1_numTips = length(E11{1, 'TipsVolumes'}{1});
E11_K1_avgTipVol = mean(E11{1, 'TipsVolumes'}{1}) * voxel_to_mm3;
E11_K1_numGen = E11{1, 'numGenerations'}{1};
E11_K1_SumClasses = E11_K1_TreeVol + E11_K1_TipsVol + E11_K1_CapVol + E11_K1_StrVol;

% K2
E11_K2_TotVol = E11{2, 'TotVol'}{1} * voxel_to_mm3;
E11_K2_TreeVol = E11{2, 'TreeVol'}{1} * voxel_to_mm3;
E11_K2_CapVol = E11{2, 'CapVol'}{1} * voxel_to_mm3;
E11_K2_StrVol = E11{2, 'StrVol'}{1} * voxel_to_mm3;
E11_K2_TipsVol = sum(E11{2, 'TipsVolumes'}{1}) * voxel_to_mm3;
E11_K2_numTips = length(E11{2, 'TipsVolumes'}{1});
E11_K2_avgTipVol = mean(E11{2, 'TipsVolumes'}{1}) * voxel_to_mm3;
E11_K2_numGen = E11{2, 'numGenerations'}{1};
E11_K2_SumClasses = E11_K2_TreeVol + E11_K2_TipsVol + E11_K2_CapVol + E11_K2_StrVol;

%% Extract data for MA40
% K1
MA40_K1_TotVol = MA40{1, 'TotVol'}{1} * voxel_to_mm3;
MA40_K1_TreeVol = MA40{1, 'TreeVol'}{1} * voxel_to_mm3;
MA40_K1_CapVol = MA40{1, 'CapVol'}{1} * voxel_to_mm3;
MA40_K1_StrVol = MA40{1, 'StrVol'}{1} * voxel_to_mm3;
MA40_K1_TipsVol = sum(MA40{1, 'TipsVolumes'}{1}) * voxel_to_mm3;
MA40_K1_numTips = length(MA40{1, 'TipsVolumes'}{1});
MA40_K1_avgTipVol = mean(MA40{1, 'TipsVolumes'}{1}) * voxel_to_mm3;
MA40_K1_numGen = MA40{1, 'numGenerations'}{1};
MA40_K1_SumClasses = MA40_K1_TreeVol + MA40_K1_TipsVol + MA40_K1_CapVol + MA40_K1_StrVol;

% K2
MA40_K2_TotVol = MA40{2, 'TotVol'}{1} * voxel_to_mm3;
MA40_K2_TreeVol = MA40{2, 'TreeVol'}{1} * voxel_to_mm3;
MA40_K2_CapVol = MA40{2, 'CapVol'}{1} * voxel_to_mm3;
MA40_K2_StrVol = MA40{2, 'StrVol'}{1} * voxel_to_mm3;
MA40_K2_TipsVol = sum(MA40{2, 'TipsVolumes'}{1}) * voxel_to_mm3;
MA40_K2_numTips = length(MA40{2, 'TipsVolumes'}{1});
MA40_K2_avgTipVol = mean(MA40{2, 'TipsVolumes'}{1}) * voxel_to_mm3;
MA40_K2_numGen = MA40{2, 'numGenerations'}{1};
MA40_K2_SumClasses = MA40_K2_TreeVol + MA40_K2_TipsVol + MA40_K2_CapVol + MA40_K2_StrVol;

%% Calculate averages
E11_avg_TotVol = mean([E11_K1_TotVol, E11_K2_TotVol]);
E11_avg_TreeVol = mean([E11_K1_TreeVol, E11_K2_TreeVol]);
E11_avg_TipsVol = mean([E11_K1_TipsVol, E11_K2_TipsVol]);
E11_avg_CapVol = mean([E11_K1_CapVol, E11_K2_CapVol]);
E11_avg_StrVol = mean([E11_K1_StrVol, E11_K2_StrVol]);
E11_avg_numTips = mean([E11_K1_numTips, E11_K2_numTips]);
E11_avg_avgTipVol = mean([E11_K1_avgTipVol, E11_K2_avgTipVol]);
E11_avg_numGen = mean([E11_K1_numGen, E11_K2_numGen]);

MA40_avg_TotVol = mean([MA40_K1_TotVol, MA40_K2_TotVol]);
MA40_avg_TreeVol = mean([MA40_K1_TreeVol, MA40_K2_TreeVol]);
MA40_avg_TipsVol = mean([MA40_K1_TipsVol, MA40_K2_TipsVol]);
MA40_avg_CapVol = mean([MA40_K1_CapVol, MA40_K2_CapVol]);
MA40_avg_StrVol = mean([MA40_K1_StrVol, MA40_K2_StrVol]);
MA40_avg_numTips = mean([MA40_K1_numTips, MA40_K2_numTips]);
MA40_avg_avgTipVol = mean([MA40_K1_avgTipVol, MA40_K2_avgTipVol]);
MA40_avg_numGen = mean([MA40_K1_numGen, MA40_K2_numGen]);

%% Calculate percentages
E11_K1_TreePct = (E11_K1_TreeVol / E11_K1_SumClasses) * 100;
E11_K1_TipsPct = (E11_K1_TipsVol / E11_K1_SumClasses) * 100;
E11_K1_CapPct = (E11_K1_CapVol / E11_K1_SumClasses) * 100;
E11_K1_StrPct = (E11_K1_StrVol / E11_K1_SumClasses) * 100;

E11_K2_TreePct = (E11_K2_TreeVol / E11_K2_SumClasses) * 100;
E11_K2_TipsPct = (E11_K2_TipsVol / E11_K2_SumClasses) * 100;
E11_K2_CapPct = (E11_K2_CapVol / E11_K2_SumClasses) * 100;
E11_K2_StrPct = (E11_K2_StrVol / E11_K2_SumClasses) * 100;

MA40_K1_TreePct = (MA40_K1_TreeVol / MA40_K1_SumClasses) * 100;
MA40_K1_TipsPct = (MA40_K1_TipsVol / MA40_K1_SumClasses) * 100;
MA40_K1_CapPct = (MA40_K1_CapVol / MA40_K1_SumClasses) * 100;
MA40_K1_StrPct = (MA40_K1_StrVol / MA40_K1_SumClasses) * 100;

MA40_K2_TreePct = (MA40_K2_TreeVol / MA40_K2_SumClasses) * 100;
MA40_K2_TipsPct = (MA40_K2_TipsVol / MA40_K2_SumClasses) * 100;
MA40_K2_CapPct = (MA40_K2_CapVol / MA40_K2_SumClasses) * 100;
MA40_K2_StrPct = (MA40_K2_StrVol / MA40_K2_SumClasses) * 100;

E11_avg_TreePct = mean([E11_K1_TreePct, E11_K2_TreePct]);
E11_avg_TipsPct = mean([E11_K1_TipsPct, E11_K2_TipsPct]);
E11_avg_CapPct = mean([E11_K1_CapPct, E11_K2_CapPct]);
E11_avg_StrPct = mean([E11_K1_StrPct, E11_K2_StrPct]);

MA40_avg_TreePct = mean([MA40_K1_TreePct, MA40_K2_TreePct]);
MA40_avg_TipsPct = mean([MA40_K1_TipsPct, MA40_K2_TipsPct]);
MA40_avg_CapPct = mean([MA40_K1_CapPct, MA40_K2_CapPct]);
MA40_avg_StrPct = mean([MA40_K1_StrPct, MA40_K2_StrPct]);

%% Create comprehensive figure with all plots
figure('Position', [100, 100, 1800, 400]);

% Define x-axis positions for individual samples with offsets
x_offset = 0.08;  % Amount of horizontal jitter
x_E11_K1 = 1 - x_offset;
x_E11_K2 = 1 + x_offset;
x_MA40_K1 = 2 - x_offset;
x_MA40_K2 = 2 + x_offset;

% Subplot 1: Overall Kidney Volume
subplot(1,5,1);
hold on;
b = bar([1, 2], [E11_avg_TotVol, MA40_avg_TotVol], 'FaceColor', 'flat', 'FaceAlpha', 0.6, 'LineWidth', 1);
b.CData = [217, 189, 213; 205, 92, 92]/255;
b.EdgeColor = 'flat';
scatter([1, 1], [E11_K1_TotVol, E11_K2_TotVol], 30, 'k', 'filled');
scatter([2, 2], [MA40_K1_TotVol, MA40_K2_TotVol], 30, 'k', 'filled');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Mouse', 'Macaque'}, 'FontSize', 12);
ylabel('Volume (mm^3)', 'FontSize', 12);
title('Metanephros Volume', 'FontSize', 12);
set(gca,'FontName', 'Arial','LineWidth', 1)
box on;
hold off;

% Subplot 2: Tissue Composition (Stacked Bar) - MODIFIED TO SHOW INDIVIDUAL SAMPLES
subplot(1,5,2);
hold on;

% Data for all 4 samples
stackedData = [E11_K1_TreePct, E11_K1_TipsPct, E11_K1_CapPct, E11_K1_StrPct;
               E11_K2_TreePct, E11_K2_TipsPct, E11_K2_CapPct, E11_K2_StrPct;
               MA40_K1_TreePct, MA40_K1_TipsPct, MA40_K1_CapPct, MA40_K1_StrPct;
               MA40_K2_TreePct, MA40_K2_TipsPct, MA40_K2_CapPct, MA40_K2_StrPct];

b = bar(stackedData, 'stacked', 'LineWidth', 1);
b(1).FaceColor = cmap_base(1,:);  % Tree
b(1).FaceAlpha = 0.6;
b(1).EdgeColor = cmap_base(1,:);
b(2).FaceColor = cmap_base(2,:);  % Vesicles
b(2).FaceAlpha = 0.6;
b(2).EdgeColor = cmap_base(2,:);
b(3).FaceColor = cmap_base(3,:);  % Cap
b(3).FaceAlpha = 0.6;
b(3).EdgeColor = cmap_base(3,:);
b(4).FaceColor = cmap_base(4,:);  % Stroma
b(4).FaceAlpha = 0.6;
b(4).EdgeColor = cmap_base(4,:);

% Add vertical line to separate species
plot([2.5 2.5], [0 100], 'k--', 'LineWidth', 0.5);

set(gca, 'XTick', [1.5 3.5], 'XTickLabel', {'Mouse', 'Macaque'}, 'FontSize', 12);
xtickangle(45);
ylabel('Volume Fraction (%)', 'FontSize', 12);
title('Tissue Composition', 'FontSize', 12);
legend('Ureteric Tree', 'Renal Vesicles', 'Cap Mesenchyme', 'Stroma', 'Location', 'best', 'FontSize', 10);
ylim([0 100]);
set(gca,'FontName', 'Arial','LineWidth', 1)
box on;
hold off;

% Subplot 3: Maximum Number of Generations - WITH X OFFSET
subplot(1,5,3);
hold on;
b = bar([1, 2], [E11_avg_numGen, MA40_avg_numGen], 'FaceColor', 'flat', 'FaceAlpha', 0.6, 'LineWidth', 1);
b.CData = [cmap_mouse(1,:); cmap_macaque(1,:)];
b.EdgeColor = 'flat';
scatter([x_E11_K1, x_E11_K2], [E11_K1_numGen, E11_K2_numGen], 30, 'k', 'filled');
scatter([x_MA40_K1, x_MA40_K2], [MA40_K1_numGen, MA40_K2_numGen], 30, 'k', 'filled');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Mouse', 'Macaque'}, 'FontSize', 12);
ylabel('Max Strahler Order', 'FontSize', 12);
title('Ureteric Tree Strahler Order', 'FontSize', 12);
set(gca,'FontName', 'Arial','LineWidth', 1)
box on;
hold off;

% Subplot 4: Number of Tips - WITH X OFFSET
subplot(1,5,4);
hold on;
b = bar([1, 2], [E11_avg_numTips, MA40_avg_numTips], 'FaceColor', 'flat', 'FaceAlpha', 0.6, 'LineWidth', 1);
b.CData = [cmap_mouse(2,:); cmap_macaque(2,:)];
b.EdgeColor = 'flat';
scatter([x_E11_K1, x_E11_K2], [E11_K1_numTips, E11_K2_numTips], 30, 'k', 'filled');
scatter([x_MA40_K1, x_MA40_K2], [MA40_K1_numTips, MA40_K2_numTips], 30, 'k', 'filled');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Mouse', 'Macaque'}, 'FontSize', 12);
ylabel('Number of Tips', 'FontSize', 12);
title('Number of Renal Vesicles', 'FontSize', 12);
set(gca,'FontName', 'Arial','LineWidth', 1)
box on;
hold off;

% Subplot 5: Average Tip Volume - WITH X OFFSET
subplot(1,5,5);
hold on;
b = bar([1, 2], [E11_avg_avgTipVol, MA40_avg_avgTipVol], 'FaceColor', 'flat', 'FaceAlpha', 0.6, 'LineWidth', 1);
b.CData = [cmap_mouse(2,:); cmap_macaque(2,:)];
b.EdgeColor = 'flat';
scatter([x_E11_K1, x_E11_K2], [E11_K1_avgTipVol, E11_K2_avgTipVol], 30, 'k', 'filled');
scatter([x_MA40_K1, x_MA40_K2], [MA40_K1_avgTipVol, MA40_K2_avgTipVol], 30, 'k', 'filled');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Mouse', 'Macaque'}, 'FontSize', 12);
ylabel('Average Tip Volume (mm^3)', 'FontSize', 12);
title('Renal Vesicle Volume', 'FontSize', 12);
set(gca,'FontName', 'Arial','LineWidth', 1)
box on;
hold off;

%% Print detailed summary statistics
fprintf('\n========== KIDNEY ANALYSIS SUMMARY ==========\n\n');

fprintf('--- MOUSE (E11) AVERAGE ---\n');
fprintf('Total Volume: %.4f mm^3\n', E11_avg_TotVol);
fprintf('Ureteric Tree: %.4f mm^3 (%.1f%%)\n', E11_avg_TreeVol, E11_avg_TreePct);
fprintf('Renal Vesicles: %.4f mm^3 (%.1f%%)\n', E11_avg_TipsVol, E11_avg_TipsPct);
fprintf('Cap Mesenchyme: %.4f mm^3 (%.1f%%)\n', E11_avg_CapVol, E11_avg_CapPct);
fprintf('Stroma: %.4f mm^3 (%.1f%%)\n', E11_avg_StrVol, E11_avg_StrPct);
fprintf('Number of Tips: %.1f\n', E11_avg_numTips);
fprintf('Average Tip Volume: %.6f mm^3\n', E11_avg_avgTipVol);
fprintf('Max Generations: %.1f\n\n', E11_avg_numGen);

fprintf('--- MACAQUE (MA40) AVERAGE ---\n');
fprintf('Total Volume: %.4f mm^3\n', MA40_avg_TotVol);
fprintf('Ureteric Tree: %.4f mm^3 (%.1f%%)\n', MA40_avg_TreeVol, MA40_avg_TreePct);
fprintf('Renal Vesicles: %.4f mm^3 (%.1f%%)\n', MA40_avg_TipsVol, MA40_avg_TipsPct);
fprintf('Cap Mesenchyme: %.4f mm^3 (%.1f%%)\n', MA40_avg_CapVol, MA40_avg_CapPct);
fprintf('Stroma: %.4f mm^3 (%.1f%%)\n', MA40_avg_StrVol, MA40_avg_StrPct);
fprintf('Number of Tips: %.1f\n', MA40_avg_numTips);
fprintf('Average Tip Volume: %.6f mm^3\n', MA40_avg_avgTipVol);
fprintf('Max Generations: %.1f\n', MA40_avg_numGen);