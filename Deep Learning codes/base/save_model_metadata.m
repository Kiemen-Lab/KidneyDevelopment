function save_model_metadata(pthDL,pthim,WS,nm,umpix,cmap,sxy,classNames,ntrain,nvalidate)
if ~isfolder(pthDL);mkdir(pthDL);end
datafile=[pthDL,'net.mat'];
disp('Saving model metadata and classification colormap...')

if strcmp(classNames(end),"black")==0
    classNames(end+1)="black";
end

if strcmp(classNames(end),"black")
    classNames=classNames(1:end-1);
end

% fix WS and classNames if there are classes to delete
ncombine=WS{3};
nload=WS{4};
ndelete=WS{5};
ndelete=sort(ndelete,'descend');
if ~isempty(ndelete)
    for b=ndelete
        oldnum=ncombine(b);
        ncombine(b)=1;
        ncombine(ncombine>oldnum)=ncombine(ncombine>oldnum)-1;
        nload=setdiff(nload,b,'stable');
        if length(classNames)==max(WS{3})
            zz=intersect(1:length(classNames),[1:oldnum-1 oldnum+1:length(classNames)]);
            classNames=classNames(zz);
            cmap=cmap(zz,:);
        end
    end
    WS{3}=ncombine;
    WS{4}=nload;
end

nwhite=WS{3};
nwhite=nwhite(WS{2});
nwhite=nwhite(1);

if max(WS{3})~=length(classNames)
    error(' The length of classNames does not match the number of classes specified in WS{3}.')
end

if strcmp(classNames(end),"black")==0
    classNames(end+1)="black";
end
nblack=length(classNames);

if exist(datafile,'file')
    variableInfo = who('-file',datafile);
    if ismember('net', variableInfo) && ismember('pthim', variableInfo)
        error([' A network has already been trained for model ',nm,...
            '. Choose a new model name to retrain.'])
    elseif ismember('net', variableInfo) && ~ismember('pthim', variableInfo)
        save([pthDL,'net.mat'],'pthim','WS','nm','umpix','cmap','sxy','nblack','nwhite','classNames','ntrain','nvalidate','-append');
        error([' A network has already been trained for model ',nm,...
            '. Metadata added to net.mat file. Choose a new model name to retrain.'])
    end

    % things that don't affect loading of xml data: nblack, sxy, cmap, classNames, nm
    save([pthDL,'net.mat'],'pthim','WS','nm','umpix','cmap','sxy','nblack','nwhite','classNames','ntrain','nvalidate','-append');
else
    % if file doesn't exist, save all the variables
    save([pthDL,'net.mat'],'pthim','WS','nm','umpix','cmap','sxy','nblack','nwhite','classNames','ntrain','nvalidate');
end


% plot color legend
plot_cmap_legend(cmap,classNames);
h=gcf;
savefig(h,[pthDL,'model_color_legend.fig']);
print(h,[pthDL,'model_color_legend'],'-djpeg','-r300');