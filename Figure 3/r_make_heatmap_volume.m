pth = '\\Lucie Dequiedt\Kidney Project\Data for paper\Figure 4\Glom Volumes\volgloms\';

for sample = {'Mac_d_labelled_glom.mat','Hum_K3_top_labelled_glom.mat','E17_K2_labelled_glom.mat'}

        load(fullfile(pth, sample{1,1}));
        
        volp = volglom;
        volglom = logical(volglom);
        
        % Get volume of each object
        V = regionprops3(volp,'Volume');
        V = cat(2, V.Volume); 
        V = V * 4 * 4 * 4;  % Convert to physical units
        
        % Remove glomeruli that are too big (likely merged objects)
        indices = find(V > 700000);
        volp(ismember(volp, indices)) = 0;
        volglom = volglom .* (volp > 0);
        
        indices = find(V == 0);
        volp(ismember(volp, indices)) = 0;
        volglom = volglom .* (volp > 0);
        
        % Recalculate volumes after removal
        V = regionprops3(volp, 'Volume');
        V = cat(2, V.Volume); 
        V = V * 4 * 4 * 4;
        
        % Sort volumes in ascending order and get the sorting indices
        [V_sorted, sort_idx] = sort(V, 'ascend');
        
        % Create new volume with labels 1->N based on volume (smallest=1, largest=N)
        voln = zeros(size(volp));
        % rearrange the cell packing density in ascending order
        %d2=[0 sort_idx.']; % cell packing densities of each glomeruli in the volume 
        d2=[0 V.'];
        voln=d2(volp+1); % matrix where glomeruli has its volume as a value  
        
        % Now voln contains glomeruli labeled 1 (smallest volume) to max (largest volume)
        % V_sorted contains the corresponding volumes in ascending order
        p1=prctile(V(V~=0),1);
        p2=prctile(V,99);
        
        voln(voln>0 & voln<p1)=p1;
        voln(voln>p2)=p2;
        save(['voln_',sample{1,1}],'voln','-v7.3')
end