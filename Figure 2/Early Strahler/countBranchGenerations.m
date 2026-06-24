% Convert skeleton to graph structure
function [maxStrahlerOrder, strahlerLabels] = countBranchGenerations(skeleton, rootPoint)
    % rootPoint: [x, y, z] coordinates of the tree root (e.g., trunk base)
    % Computes Strahler order for each point in the skeleton
    
    % Get all skeleton points
    [skelX, skelY, skelZ] = ind2sub(size(skeleton), find(skeleton));
    skelPoints = [skelX, skelY, skelZ];
    numPoints = size(skelPoints, 1);
    
    % Build adjacency matrix (26-connectivity)
    adjacency = sparse(numPoints, numPoints);
    
    for i = 1:numPoints
        % Find neighbors within distance sqrt(3) (26-connected)
        dists = sqrt(sum((skelPoints - skelPoints(i,:)).^2, 2));
        neighbors = find(dists > 0 & dists <= sqrt(3));
        adjacency(i, neighbors) = 1;
    end
    
    % Find root index
    rootDist = sqrt(sum((skelPoints - rootPoint).^2, 2));
    [~, rootIdx] = min(rootDist);
    
    % Identify node types
    degree = full(sum(adjacency, 2));
    endPoints = find(degree == 1);
    branchPoints = find(degree > 2);
    
    fprintf('End points: %d\n', length(endPoints));
    fprintf('Branch points: %d\n', length(branchPoints));
    
    % Initialize Strahler orders
    strahlerLabels = zeros(numPoints, 1);
    visited = false(numPoints, 1);
    
    % Assign Strahler order 1 to all endpoints
    strahlerLabels(endPoints) = 1;
    visited(endPoints) = true;
    
    % Build tree structure from root
    parent = -1 * ones(numPoints, 1);
    children = cell(numPoints, 1);
    
    % BFS to establish parent-child relationships
    queue = rootIdx;
    visited_bfs = false(numPoints, 1);
    visited_bfs(rootIdx) = true;
    
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];
        
        neighbors = find(adjacency(current, :));
        for j = 1:length(neighbors)
            neighbor = neighbors(j);
            if ~visited_bfs(neighbor)
                visited_bfs(neighbor) = true;
                parent(neighbor) = current;
                children{current} = [children{current}, neighbor];
                queue = [queue, neighbor];
            end
        end
    end
    
    % Compute Strahler order recursively from leaves to root
    strahlerLabels = computeStrahlerRecursive(rootIdx, children, degree, strahlerLabels);
    
    maxStrahlerOrder = max(strahlerLabels);
    
    % Visualize Strahler orders
    figure;
    scatter3(skelPoints(:,1), skelPoints(:,2), skelPoints(:,3), ...
        20, strahlerLabels, 'filled');
    colorbar;
    colormap(jet);
    title(sprintf('Strahler Order (Max: %d)', maxStrahlerOrder));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;
    
    % Display statistics
    fprintf('\nStrahler Order Distribution:\n');
    for order = 1:maxStrahlerOrder
        count = sum(strahlerLabels == order);
        fprintf('  Order %d: %d points (%.1f%%)\n', order, count, 100*count/numPoints);
    end
end

function strahlerOrder = computeStrahlerRecursive(nodeIdx, children, degree, strahlerOrder)
    % Recursive function to compute Strahler order
    
    % If this is an endpoint, it already has order 1
    if degree(nodeIdx) == 1 && isempty(children{nodeIdx})
        strahlerOrder(nodeIdx) = 1;
        return;
    end
    
    % Get children
    childList = children{nodeIdx};
    
    if isempty(childList)
        % Leaf node
        strahlerOrder(nodeIdx) = 1;
        return;
    end
    
    % Recursively compute Strahler order for all children
    childOrders = zeros(length(childList), 1);
    for i = 1:length(childList)
        strahlerOrder = computeStrahlerRecursive(childList(i), children, degree, strahlerOrder);
        childOrders(i) = strahlerOrder(childList(i));
    end
    
    % Apply Strahler ordering rules:
    % 1. If all children have different orders, take the maximum
    % 2. If at least two children have the same maximum order, increment by 1
    
    maxOrder = max(childOrders);
    numMaxOrder = sum(childOrders == maxOrder);
    
    if numMaxOrder >= 2
        % Two or more children with same max order
        strahlerOrder(nodeIdx) = maxOrder + 1;
    else
        % All children have different orders
        strahlerOrder(nodeIdx) = maxOrder;
    end
end