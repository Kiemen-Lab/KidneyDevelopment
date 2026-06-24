pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};
noise_label = 14;

classNames = {'Urothelium','Med. Coll. Duct','Coll. Duct','Dist. Tub.','Henle Loop','Prox. Tub.', ...
                  'Glom. Tuft','Bow. Cap.','Vein','Arteries','Arterioles','Stroma','Dev. Corpuscle','Dev. Nephron','Blastema'};
%% Compute Cell Fractions
cellfrac = zeros(8,16);  % Changed from 10 to 8 rows
load('diameters_nuclei.mat')
sz=4;

% Temporary storage for cell counts before normalization
cellcounts = zeros(8,16);

for i=1:length(samples)
    load([pth,samples{i},'.mat'],'volTA');volTA=double(volTA);
    load([pth,samples{i},'_cells.mat'],'vc');vc = double(vc);
    volTA(volTA==noise_label) = 0;
    
    temp_counts = zeros(1,16);
    for j=1:length(classNames)+1
        cells = vc.*(volTA==j);
        cells = sum(cells(:));
        temp_counts(j) = cells* (sz/(sz+D(j)));
    end
    
    % Combine macaque samples 2 by 2
    if i == 8 || i == 9  % Mac_b and Mac_c combine into row 7
        cellcounts(i-1,:) = cellcounts(i-1,:) + temp_counts;
    elseif i == 10  % Mac_d combines into row 8
        cellcounts(i-2,:) = cellcounts(i-2,:) + temp_counts;
    else  % E17_K1, E17_K2, Hum_K1, Hum_K2, Hum_K3_bottom, Hum_K3_top, Mac_a
        cellcounts(i,:) = temp_counts;
    end
end

% Normalize after combining
for i=1:8
    cellfrac(i,:) = (cellcounts(i,:)/sum(cellcounts(i,:)))*100;
end

save([pth,'cell_fraction.mat'],"cellfrac");