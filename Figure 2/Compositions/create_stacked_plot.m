function create_stacked_plot(compo,cmap,classNames,sampNames,name_species,displeg)

    fontSize = 14;
    numSamples  = size(compo,1);
    numElements = size(compo,2);
    
    figure;
    bar(compo, 'stacked');
    for i = 1:numElements
        num = numElements-i+1;
        b = findobj(gca,'Type','Bar');
        b(i).FaceColor = cmap(num,:);
        b(i).FaceAlpha = 0.7;
        b(i).EdgeColor = cmap(num,:);
    end
    % Axis formatting
    ylabel('Composition (%)');
    xticks(1:numSamples);
    xticklabels(sampNames);
    ylim([0 100]);
    
    % Legend
    if displeg == 1
        legend(classNames, 'Location', 'eastoutside');
    end
    
    set(gca,'FontSize',12);
    box off;

    set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', fontSize, ...
    'Box', 'on', ...
    'LineWidth', 1);

    ylabel('Composition [%]', ...
    'FontName', 'Arial', ...
    'FontSize', fontSize);

    title(name_species, ...
    'FontName', 'Arial', ...
    'FontSize', fontSize, ...
    'FontWeight', 'normal');

end