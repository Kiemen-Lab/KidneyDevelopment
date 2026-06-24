path(path,'\\lucie\codes_LD\Fetal Kidney\figures\z projections\')
cmap = [191   188   109; ...       % 5 distal tubule
        135   214   193; ...       % 6 proximal tubule
        113   191   109; ...       % 7 henle Loop
        235   186   134];      % 8 collecting duct
titles = {"Distal Tubule" "Proximal Tubule" "Loop of Henle" "Collecting Duct"};
make_cmap_legend(cmap,titles);

%%
cmap = [224, 168, 194; ...       % 5 distal tubule
        235   154   134; ...       % 6 proximal tubule
        235   186   134];      % 8 collecting duct
titles = {"Urothelium" "Medullary Collecting Duct" "Collecting Duct"};
make_cmap_legend(cmap,titles);