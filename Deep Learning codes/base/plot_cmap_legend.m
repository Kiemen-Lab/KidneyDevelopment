function plot_cmap_legend(cmap,titles)

im=uint8([]);
for k=1:size(cmap,1)
    tmp=permute(cmap(k,:),[1 3 2]);
    tmp=repmat(tmp,[50 50 1]);
    im=cat(2,im,tmp);
end

if size(cmap,1)==length(titles)-1
    titles=titles(1:end-1);
end

for b=1:length(titles);titles(b)=strrep(titles(b),'_',' ');end

if exist('titles','var')
    im=imrotate(im,270);
    figure(182);imagesc(im)
    axis equal
    xlim([0 size(im,2)]);ylim([0 size(im,1)])
    yticks(15:50:size(im,1))
    set(gca,'TickLength',[0 0])
    xticks([])
    yticklabels(titles)%ytickangle(90)
    set(gca,'fontsize',15);
else
    figure(182);imshow(im)
end

set(gcf,'color','w');set(gca,'LineWidth',1);
pause(0.05)
