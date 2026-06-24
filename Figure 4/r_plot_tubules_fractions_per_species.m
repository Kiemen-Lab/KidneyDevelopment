load("\\Lucie Dequiedt\Kidney Project\Volumes\tubule_analysis\tubule_analysis_results.mat")

config.volumeThresholds = struct(...
    'E17_K1', 0.1, 'E17_K2', 0.1, ...
    'Mac_a', 0.4, 'Mac_b', 0.4, ...
    'Hum_K1', 0.2, 'Hum_K2', 0.2, ...
    'Hum_K3_bottom', 0.4, 'Hum_K3_top', 0.4);

cmap = [191   188   109; ...       % 5 distal tubule
        135   214   193; ...       % 6 proximal tubule
        113   191   109; ...       % 7 henle Loop
        235   186   134]/255;      % 8 collecting duct

tubule_names = {'Distal tubule', 'Proximal tubule', 'Henle Loop', 'Collecting duct'};

% Create figure with 3 subplots
figure('Position', [100 100 1400 400]);

%% Mouse (subplot 1)
subplot(3,1,1)
% Get data for mouse samples
dis_mouse(:,:,1) = results.E17_K1.tubuleFractions*100;
dis_mouse(:,:,2) = results.E17_K2.tubuleFractions*100;

% Get volumes and bin centers
vols_mouse{1} = results.E17_K1.volumes * 4*4*4*10^-9;
vols_mouse{2} = results.E17_K2.volumes * 4*4*4*10^-9;
bin_mouse = results.E17_K1.binCenters * 4/1000;

% Apply volume thresholding and plot
for i = 1:4
    % Sample 1
    valid_bins = vols_mouse{1} > config.volumeThresholds.E17_K1;
    plot(bin_mouse(valid_bins), dis_mouse(valid_bins,i,1), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off','MarkerFaceColor',cmap(i,:));
    hold on
    
    % Sample 2
    valid_bins = vols_mouse{2} > config.volumeThresholds.E17_K2;
    plot(bin_mouse(valid_bins), dis_mouse(valid_bins,i,2), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'DisplayName', tubule_names{i},'MarkerFaceColor',cmap(i,:));
end
xlim([0 1])
title('Mouse')
legend('Location', 'best')
hold off

%% Human (subplot 2)
subplot(3,1,2)
% Get data for human samples
dis_human(:,:,1) = results.Hum_K1.tubuleFractions*100;
dis_human(:,:,2) = results.Hum_K2.tubuleFractions*100;
dis_human(:,:,3) = results.Hum_K3_top.tubuleFractions*100;
dis_human(:,:,4) = results.Hum_K3_bottom.tubuleFractions*100;

% Get volumes and bin centers
vols_human{1} = results.Hum_K1.volumes * 4*4*4*10^-9;
vols_human{2} = results.Hum_K2.volumes * 4*4*4*10^-9;
vols_human{3} = results.Hum_K3_top.volumes * 4*4*4*10^-9;
vols_human{4} = results.Hum_K3_bottom.volumes * 4*4*4*10^-9;
bin_human = results.Hum_K1.binCenters * 4/1000;

% Apply volume thresholding and plot
for i = 1:4
    % Sample 1
    valid_bins = vols_human{1} > config.volumeThresholds.Hum_K1;
    plot(bin_human(valid_bins), dis_human(valid_bins,i,1), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off','MarkerFaceColor',cmap(i,:));
    hold on
    
    % Sample 2
    valid_bins = vols_human{2} > config.volumeThresholds.Hum_K2;
    plot(bin_human(valid_bins), dis_human(valid_bins,i,2), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off','MarkerFaceColor',cmap(i,:));
    
    % Sample 3
    valid_bins = vols_human{3} > config.volumeThresholds.Hum_K3_top;
    plot(bin_human(valid_bins), dis_human(valid_bins,i,3), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off','MarkerFaceColor',cmap(i,:));
    
    % Sample 4
    valid_bins = vols_human{4} > config.volumeThresholds.Hum_K3_bottom;
    plot(bin_human(valid_bins), dis_human(valid_bins,i,4), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off');
end
xlim([0 1])
ylabel('Volume Fraction [%]')
title('Human')
hold off

%% Macaque (subplot 3)
subplot(3,1,3)
% Combine Mac_a and Mac_b
counts_mac1 = results.Mac_a.tubuleCounts + results.Mac_b.tubuleCounts;
vols_mac1 = results.Mac_a.volumes + results.Mac_b.volumes;
dis_mac(:,:,1) = counts_mac1./vols_mac1*100;

% Combine Mac_c and Mac_d
counts_mac2 = results.Mac_c.tubuleCounts + results.Mac_d.tubuleCounts;
vols_mac2 = results.Mac_c.volumes + results.Mac_d.volumes;
dis_mac(:,:,2) = counts_mac2./vols_mac2*100;

% Get volumes in mm^3 for thresholding
vols_mac_thresh{1} = vols_mac1 * 4*4*4*10^-9;
vols_mac_thresh{2} = vols_mac2 * 4*4*4*10^-9;
bin_mac = results.Mac_a.binCenters * 4/1000;

% Apply volume thresholding and plot
for i = 1:4
    % Sample 1 (Mac_a + Mac_b)
    valid_bins = vols_mac_thresh{1} > config.volumeThresholds.Mac_a;
    plot(bin_mac(valid_bins), dis_mac(valid_bins,i,1), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off','MarkerFaceColor',cmap(i,:));
    hold on
    
    % Sample 2 (Mac_c + Mac_d)
    valid_bins = vols_mac_thresh{2} > config.volumeThresholds.Mac_b;
    plot(bin_mac(valid_bins), dis_mac(valid_bins,i,2), '-o', 'Color', cmap(i,:), 'LineWidth', 2, 'HandleVisibility', 'off','MarkerFaceColor',cmap(i,:));
end
xlim([0 1])
title('Macaque')
hold off

% Add common x-axis label
han = axes(gcf, 'visible', 'off');
han.XLabel.Visible = 'on';
xlabel(han, 'Distance to kidney surface [mm]');