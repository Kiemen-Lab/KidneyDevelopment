function xyout=load_xml_coords(xmlfile)

str = fileread(xmlfile);

% extract starting text
bb=strfind(str,'<Annotation Id=');
toptext=str(1:bb(1)-1);
str=str(bb:end);

bb=strfind(str,'<Annotation Id=');
annlayers=cell([1 length(bb)]);
for k=1:length(bb)
    num=bb(k);
    if k==length(bb)
        txt=str(num:end);
        Annfind=strfind(txt,'</Annotations>');
        bottext=txt(Annfind:end);
        txt=txt(1:Annfind-1);
        annlayers{k}=txt;
    else
        num2=bb(k+1);
        annlayers{k}=str(num:num2-1);
    end
end

% extract coordinates from annlayers
xyout=[];
for b=1:length(annlayers)
    txt=annlayers{b};
    numReg=strfind(txt,'<Region Id=');
    numReg=cat(2,numReg,length(txt));
    % isolate text for each region
    for z=1:length(numReg)-1
        num=numReg([z z+1]);
        txtReg=txt(num(1):num(2));
        fndxyA=strfind(txtReg,'Vertex X=');
        fndxy=cat(2,fndxy,length(txtReg));
        %fndxyB=strfind(txtReg,' Y=');
        %fndxyC=strfind(txtReg,' Z=');
        for pts=1:length(fndxyA)-1
            numpt=[fndxyA(pts) fndxyB(pts) fndxyC(pts)];
            lineXY=txtReg(fndxy(pts):fndxy(pts+1));
            %lineX=str2double(txtReg(numpt(1)+10:numpt(2)-2));
            %lineY=str2double(txtReg(numpt(2)+4:numpt(3)-2));

            lineX=[];lineY=[];
            vals=[b z lineX lineY];
            xyout=cat(1,xyout,vals);
        end
        disp([z length(numReg)])
    end
end



end