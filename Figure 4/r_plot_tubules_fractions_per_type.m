load("\\Lucie Dequiedt\Kidney Project\Volumes\tubule_analysis\tubule_analysis_results.mat")

config.volumeThresholds = struct(...
    'E17_K1', 0.1, 'E17_K2', 0.1, ...
    'Mac_a', 0.4, 'Mac_b', 0.4, ...
    'Hum_K1', 0.2, 'Hum_K2', 0.2, ...
    'Hum_K3_bottom', 0.4, 'Hum_K3_top', 0.4);

% Tubule colors
cmap = [191   188   109; ...       % 5 distal tubule
        135   214   193; ...       % 6 proximal tubule
        113   191   109; ...       % 7 henle Loop
        235   186   134]/255;      % 8 collecting duct

tubule_names = {'Distal tubule', 'Proximal tubule', 'Loop of Henle', 'Collecting duct'};

% Species definition and colors
species(1).name = 'Mouse';
species(1).pattern = 'E17';

species(2).name = 'Human';
species(2).pattern = 'Hum';

species(3).name = 'Macaque';
species(3).pattern = 'Mac';

cols = [217, 189, 213; ...
        130, 73, 73; ...
        205, 92, 92]/255;

%% Prepare Mouse data
dis_mouse(:,:,1) = results.E17_K1.tubuleFractions*100;
dis_mouse(:,:,2) = results.E17_K2.tubuleFractions*100;
vols_mouse{1} = results.E17_K1.volumes * 4*4*4*10^-9;
vols_mouse{2} = results.E17_K2.volumes * 4*4*4*10^-9;
bin_mouse = results.E17_K1.binCenters * 4/1000;

%% Prepare Human data
dis_human(:,:,1) = results.Hum_K1.tubuleFractions*100;
dis_human(:,:,2) = results.Hum_K2.tubuleFractions*100;
dis_human(:,:,3) = results.Hum_K3_top.tubuleFractions*100;
dis_human(:,:,4) = results.Hum_K3_bottom.tubuleFractions*100;
vols_human{1} = results.Hum_K1.volumes * 4*4*4*10^-9;
vols_human{2} = results.Hum_K2.volumes * 4*4*4*10^-9;
vols_human{3} = results.Hum_K3_top.volumes * 4*4*4*10^-9;
vols_human{4} = results.Hum_K3_bottom.volumes * 4*4*4*10^-9;
bin_human = results.Hum_K1.binCenters * 4/1000;

%% Prepare Macaque data
counts_mac1 = results.Mac_a.tubuleCounts + results.Mac_b.tubuleCounts;
vols_mac1 = results.Mac_a.volumes + results.Mac_b.volumes;
dis_mac(:,:,1) = counts_mac1./vols_mac1*100;

counts_mac2 = results.Mac_c.tubuleCounts + results.Mac_d.tubuleCounts;
vols_mac2 = results.Mac_c.volumes + results.Mac_d.volumes;
dis_mac(:,:,2) = counts_mac2./vols_mac2*100;

vols_mac_thresh{1} = vols_mac1 * 4*4*4*10^-9;
vols_mac_thresh{2} = vols_mac2 * 4*4*4*10^-9;
bin_mac = results.Mac_a.binCenters * 4/1000;

%% Create figure with 4 subplots (one per tubule type)
figure('Position', [100 100 1400 1000]);

for tubule_idx = 1:4
    subplot(4, 1, tubule_idx)
    hold on
    
    %% Plot Macaque data FIRST
    for sample_idx = 1:2
        valid_bins = vols_mac_thresh{sample_idx} > config.volumeThresholds.(['Mac_' char('a' + sample_idx - 1)]);
        if sample_idx == 1
            plot(bin_mac(valid_bins), dis_mac(valid_bins, tubule_idx, sample_idx), ...
                '-o', 'Color', cols(3,:), 'LineWidth', 2, ...
                'MarkerFaceColor', cols(3,:), 'DisplayName', species(3).name,'MarkerSize',2);
        else
            plot(bin_mac(valid_bins), dis_mac(valid_bins, tubule_idx, sample_idx), ...
                '-o', 'Color', cols(3,:), 'LineWidth', 2, ...
                'MarkerFaceColor', cols(3,:), 'HandleVisibility', 'off','MarkerSize',2);
        end
    end
    
    %% Plot Human data SECOND
    human_samples = {'Hum_K1', 'Hum_K2', 'Hum_K3_top', 'Hum_K3_bottom'};
    for sample_idx = 1:4
        valid_bins = vols_human{sample_idx} > config.volumeThresholds.(human_samples{sample_idx});
        if sample_idx == 1
            plot(bin_human(valid_bins), dis_human(valid_bins, tubule_idx, sample_idx), ...
                '-o', 'Color', cols(2,:), 'LineWidth', 2, ...
                'MarkerFaceColor', cols(2,:), 'DisplayName', species(2).name,'MarkerSize',2);
        else
            plot(bin_human(valid_bins), dis_human(valid_bins, tubule_idx, sample_idx), ...
                '-o', 'Color', cols(2,:), 'LineWidth', 2, ...
                'MarkerFaceColor', cols(2,:), 'HandleVisibility', 'off','MarkerSize',2);
        end
    end
    
    %% Plot Mouse data LAST
    for sample_idx = 1:2
        valid_bins = vols_mouse{sample_idx} > config.volumeThresholds.(['E17_K' num2str(sample_idx)]);
        if sample_idx == 1
            plot(bin_mouse(valid_bins), dis_mouse(valid_bins, tubule_idx, sample_idx), ...
                '-o', 'Color', cols(1,:), 'LineWidth', 2, ...
                'MarkerFaceColor', cols(1,:), 'DisplayName', species(1).name,'MarkerSize',2);
        else
            plot(bin_mouse(valid_bins), dis_mouse(valid_bins, tubule_idx, sample_idx), ...
                '-o', 'Color', cols(1,:), 'LineWidth', 2, ...
                'MarkerFaceColor', cols(1,:), 'HandleVisibility', 'off','MarkerSize',2);
        end
    end
    
    %% Format subplot
    xlim([0 1])
    ylim([0 50])  % Changed to 50%
    title(tubule_names{tubule_idx}, 'FontSize', 12, 'FontWeight', 'bold')
    ylabel('Volume Fraction [%]')
    
    if tubule_idx==4
    legend('Location', 'best')
    xlabel('Distance to nephrogenic zone [mm]')
    end
    box on
    hold off
end

% Add overall title
%sgtitle('Tubule Volume Fractions Across Species', 'FontSize', 14, 'FontWeight', 'bold');