function reload_xml=check_if_model_parameters_changed(datafile,WS,umpix,nwhite,pthim)
    reload_xml=0;
    
    if ~exist(datafile,'file')
        return
    end

    % check if variables exist in datafile
    variableInfo = who('-file',datafile);
    if ~ismember('WS', variableInfo)
        disp('    WS, umpix, nwhite, and pthim do not exist in annotations.mat. Reload the xml file to add them.')
        reload_xml=1;
        return
    end

    % check if any variables have changed
    checkvar=load(datafile,'WS');checkvar=checkvar.WS;
    if ~isequal(WS,checkvar)
        disp('   reload annotation data with updated WS')
        reload_xml=1;
    end
    
    checkvar=load(datafile,'umpix');checkvar=checkvar.umpix;
    if ~isequal(umpix,checkvar)
        disp('   reload annotation data with updated umpix')
        reload_xml=1;
    end
    
    checkvar=load(datafile,'nwhite');checkvar=checkvar.nwhite;
    if ~isequal(nwhite,checkvar)
        disp('   reload annotation data with updated nwhite')
        reload_xml=1;
    end

    % checkvar=load(datafile,'pthim');checkvar=checkvar.pthim;
    % if ~isequal(pthim,checkvar)
    %     disp('   reload annotation data with updated pthim')
    %     reload_xml=1;
    % end

end