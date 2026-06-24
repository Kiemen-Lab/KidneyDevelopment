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
load([pth,'cell_fraction.mat'],"cellfrac");
%% Reorder Compositions so that the plot looks nice
% new order classes
new_order = [1,15,8,5,7,6,3,4,2,13,16,9,10,12,11,14];
compoorder = cellfrac(:,new_order);
% remove column corresponding to new position of noise value
compoorder(:,end) = [];


%% Plot Compositions Mouse 
sampNames = {'Sample A' 'Sample B'};
compo = compoorder(1:2,:);
create_stacked_plot(compo,cmap,classNames,sampNames,'Mouse',0)
ylabel('Cell Fraction [%]');
%% Plot Compositions Human 
sampNames = {'Sample A' 'Sample B' 'Sample C' 'Sample D'};
compo = compoorder(3:6,:);
create_stacked_plot(compo,cmap,classNames,sampNames,'Human',0)
ylabel('Cell Fraction [%]');
%% Plot Compositions Macaque
sampNames = {'Sample A' 'Sample B'};
compo = compoorder(7:8,:);
create_stacked_plot(compo,cmap,classNames,sampNames,'Macaque',1)
ylabel('Cell Fraction [%]');
%% Make final plot
sampNames = {'E17 (A)' 'E17 (B)' '' '' '15 weeks' '17 weeks' '19 weeks (A)' '19 weeks (B)' '' '' 'E80 (A)' 'E80 (B)'};
compo = [compoorder(1:2,:);nan(2,15);compoorder(3:6,:);nan(2,15);compoorder(7:8,:)];
create_stacked_plot(compo,cmap,classNames,sampNames,'',1)
ylabel('Cell Fraction [%]');