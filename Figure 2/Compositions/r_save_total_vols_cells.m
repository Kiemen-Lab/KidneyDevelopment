pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};
noise_label = 14;

classNames = {'Urothelium','Med. Coll. Duct','Coll. Duct','Dist. Tub.','Henle Loop','Prox. Tub.', ...
              'Glom. Tuft','Bow. Cap.','Vein','Arteries','Arterioles','Stroma','Dev. Corpuscle','Dev. Nephron','Blastema'};

load('diameters_nuclei.mat')
sz = 4;

%% Initialize storage (8 rows for combined samples)
total_volume = zeros(8, 1);      % Total volume per sample
total_cells = zeros(8, 1);       % Total cell count per sample
volume_by_class = zeros(8, 16);  % Volume per class per sample
cells_by_class = zeros(8, 16);   % Cell count per class per sample

%% Compute for each sample
for i = 1:length(samples)
    load([pth, samples{i}, '.mat'], 'volTA');
    volTA = double(volTA);
    load([pth, samples{i}, '_cells.mat'], 'vc');
    vc = double(vc);
    
    volTA(volTA == noise_label) = 0;
    
    % Determine which row to store results
    if i == 8 || i == 9  
        row_idx = i-1;
    elseif i == 10      
        row_idx = i-2;
    else                 
        row_idx = i;
    end
    
    % Volume computation
    temp_volume = histcounts(volTA, 1:17);
    volume_by_class(row_idx, :) = volume_by_class(row_idx, :) + temp_volume;
    total_volume(row_idx) = total_volume(row_idx) + sum(temp_volume);
    
    % Cell computation
    for j = 1:length(classNames)+1
        cells = vc .* (volTA == j);
        cells = sum(cells(:));
        cell_count = cells * (sz / (sz + D(j)));
        cells_by_class(row_idx, j) = cells_by_class(row_idx, j) + cell_count;
        total_cells(row_idx) = total_cells(row_idx) + cell_count;
    end
end

%% Display results
sample_labels = {'E17_K1', 'E17_K2', 'Hum_K1', 'Hum_K2', 'Hum_K3_bottom', 'Hum_K3_top', 'Mac_a+b', 'Mac_c+d'};

fprintf('\n=== TOTAL VOLUME AND CELL COUNT PER SAMPLE ===\n');
for i = 1:8
    fprintf('%s: Volume = %.9f mm^3, Cells = %.0f\n', sample_labels{i}, total_volume(i)*4*4*4*10^-9, total_cells(i));
end
%% Display results
sample_labels = {'E17_K1', 'E17_K2', 'Hum_K1', 'Hum_K2', 'Hum_K3_bottom', 'Hum_K3_top', 'Mac_a+b', 'Mac_c+d'};

fprintf('\n=== CELL DENSITY PER SAMPLE ===\n');
for i = 1:8
    fprintf('%s: Cell Density = %.2f cells/mm^3\n', sample_labels{i},total_cells(i)/(total_volume(i)*4*4*4*10^-9));
end

%% Save results
save([pth, 'volume_and_cell_summary.mat'], 'total_volume', 'total_cells', 'volume_by_class', 'cells_by_class', 'sample_labels');