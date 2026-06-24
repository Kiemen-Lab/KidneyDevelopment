path(path,'\lucie\codes_LD\Fetal Kidney\figures\z projections\')
cmap = [246 126 125;255 215 155;97 25 63];
titles = {"Bladder" "Metanephros" "Mesonephros"};
make_cmap_legend(cmap,titles);

cmap = [255 215 155;135 76 57;193 174 167];
titles = {"Ureteric Tree" "Renal Vesicle" "Metanephric Mesenchyme"};
make_cmap_legend(cmap,titles);

cmap = [97 25 63;200 121 140;242 206 189;178 130 183];
titles = {"Mesonephric Duct" "Caudal Tubules" "Primitive Glomerular Tufts" "Paramesonephric Duct"};
make_cmap_legend(cmap,titles);


%%
cmap = [219 181 194;161 194 209];
titles = {"Glomerular Tuft" "Bowmans Capsule"};
make_cmap_legend(cmap,titles);
%%
cmap = [97 25 63;246 235 250];
titles = {"Mesonephric Duct" "ECM"};
make_cmap_legend(cmap,titles);