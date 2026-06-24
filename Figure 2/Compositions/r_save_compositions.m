pth = '\\Lucie Dequiedt\Kidney Project\Volumes\';
samples = {'E17_K1' 'E17_K2' 'Hum_K1' 'Hum_K2' 'Hum_K3_bottom' 'Hum_K3_top' 'Mac_a' 'Mac_b' 'Mac_c' 'Mac_d'};
noise_label = 14;

cmap = [235   134   181; ...       % 1 ureter
        235   154   134; ...       % 15 capillary ducts
        235   186   134; ...       % 8 collecting duct
        191   188   109; ...       % 5 distal tubule
        113   191   109; ...       % 7 henle Loop
        135   214   193; ...       % 6 proximal tubule
        134   166   235; ...       % 3 glomerular tuft
        41    70   133; ...        % 4 bowman's capsule
        43    76   207; ...        % 2 vein
        145    29    29; ...       % 13 arteries
        64     3     3; ...        % 16 arterioles
        246   232   250; ...       % 9 stroma
        168   134   235; ...       % 10 developing corpuscle
        100    66   168; ...       % 12 developing nephron
        144    55   148]/255;          % 11 undifferentiated blastema cells

classNames = {'Urothelium','Med. Coll. Duct','Coll. Duct','Dist. Tub.','Henle Loop','Prox. Tub.', ...
                  'Glom. Tuft','Bow. Cap.','Vein','Arteries','Arterioles','Stroma','Dev. Corpuscle','Dev. Nephron','Blastema'};
%% Compute Compositions
compositions = zeros(8,16);

for i=1:length(samples)
    load([pth,samples{i},'.mat'],'volTA');
    load(['\\Lucie Dequiedt\Kidney Project\Volumes\inside\',samples{i},'.mat'],'volsurface');
    volTA = double(volTA).*double(volsurface);
    volTA(volTA==noise_label) = 0;
    if i == 8 | i==9
       compositions(i-1,:) = compositions(i-1,:) + histcounts(volTA,1:17);
    elseif i==10
        compositions(i-2,:) = compositions(i-2,:) + histcounts(volTA,1:17);
    else
    compositions(i,:) = histcounts(volTA,1:17);
    end
end
compositions(:,:) = (compositions(:,:)./sum(compositions(:,:),2))*100;
save(['\\Lucie Dequiedt\Kidney Project\Volumes\','composition.mat'],"compositions");
