function [numGenerations, generationLabels,skeleton] = compute_Strahler(tree,label,clean,rootPoint)
tree = tree==label;
tree = bwareaopen(tree,100);
if ~exist('clean','var')
tree = imclose(tree,strel('sphere',5));
end
tree = imfill(tree,'holes');
if exist('clean','var')
    tree = imclose(tree,strel('sphere',2));tree = imopen(tree,strel('sphere',2));
    tree = bwareaopen(tree,10000);
end
t = zeros(size(tree));
for i=1:size(tree,3)
    t(:,:,i) = imfill(tree(:,:,i),'holes');
end
tree = t;
clearvars t 

%% binaryTree = your 3D binary volume (logical array)
binaryTree = logical(tree);
% Method 1: Using built-in skeleton3D (MATLAB R2018b+)
%skeleton = bwskel(binaryTree, 'MinBranchLength', 5);
skeleton = bwskel(binaryTree);
% Method 2: Alternative using thinning (older MATLAB versions)
% skeleton = bwmorph3(binaryTree, 'skel', Inf);

% Visualize the skeleton
figure;
[x,y,z] = ind2sub(size(skeleton), find(skeleton));
plot3(x, y, z, '.', 'MarkerSize',10);
axis equal; grid on;
title('Skeleton');

%%
% Find branch points (junctions) and endpoints
branchPoints = bwmorph3(skeleton, 'branchpoints');
endPoints = bwmorph3(skeleton, 'endpoints');

% Get coordinates
[bpX, bpY, bpZ] = ind2sub(size(skeleton), find(branchPoints));
[epX, epY, epZ] = ind2sub(size(skeleton), find(endPoints));

fprintf('Branch points: %d\n', length(bpX));
fprintf('End points: %d\n', length(epX));

% Visualize
figure;
[x,y,z] = ind2sub(size(skeleton), find(skeleton));
plot3(x, y, z, 'b.', 'MarkerSize', 1); hold on;
plot3(bpX, bpY, bpZ, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot3(epX, epY, epZ, 'g^', 'MarkerSize', 10, 'LineWidth', 2);
legend('Skeleton', 'Branch Points', 'End Points');
axis equal; grid on;

ip = find(epZ==max(epZ));
if ~exist('rootPoint','var')
    rootPoint = [epX(ip) epY(ip) epZ(ip)];
end

[numGenerations, generationLabels] = countBranchGenerations(skeleton, rootPoint);
end
