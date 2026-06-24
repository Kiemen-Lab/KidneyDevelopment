pth1 = '\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 4\Glom Volumes\volgloms\';
pth2 = '\\Lucie Dequiedt\Kidney Project\Volumes\';
for sample = {'Mac_d','Hum_K3_top','E17_K2'} 

load(fullfile(pth1, [sample{1,1},'_labelled_glom.mat']));
    load(fullfile(pth2, [sample{1,1},'_cells.mat']));

    volp = volglom;
    volcell = double(vc);
    load("diameters_nuclei.mat", 'D');
    sz=4;

    volc = volp.*volcell;
    
    % for each glomeruli in the volume, compute the cell packing density 
    nummet=max(volc(:));
    distrib = histcounts(volc,0:nummet);
    distrib = distrib* (sz/(sz+D(3)));
    
    vol = histcounts(volp,0:nummet);
    vol = vol*4*4*4/10^9;
    distrib = distrib./vol;
    
    indices = find(vol>(700000/10^9)); % removing glomeruli that are too big (they are often double gloms too close together)
    % indices = find(distrib>)
    volp(ismember(volp, indices)) = 0;
    volglom=volglom.*(volp>0);
    
    % rearrange the cell packing density in ascending order
    d2=[0 distrib]; % cell packing densities of each glomeruli in the volume 
    voln=d2(volp+1); % matrix where glomeruli have their cell packing density as value 
    % what this operation does: object in matrix volp with value ex 1 will get
    % value in position 1
    
    tmp=distrib;
    p1=prctile(distrib,1);
    p2=prctile(distrib,99);
    voln(voln>p2)=p2;
    voln(voln>0 & voln<p1)=p1;

    
    voln(voln>0 & voln<p1)=p1;
    voln(voln>p2)=p2;
    save(['voln_density_',sample{1,1},'.mat'],'voln','-v7.3')

end