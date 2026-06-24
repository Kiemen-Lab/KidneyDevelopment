function [distances, volumes, cellcount, indices] = get_dist_centroid(volglom, distK, volcell, D, class)
% GET_DIST_CENTROID
% Computes the distance to the kidney edge at the centroid of each glomerulus,
% its volume, and (optionally) the count of distinct cells per glomerulus.
%
% INPUTS
%   volglom : 3D labeled volume, each glomerulus has a unique integer label
%   distK   : 3D matrix, distance of each voxel to the kidney outer edge
%   volcell : 3D binary volume, each cell ==1
%   D       : vector of mean diameters for each cell class (optional)
%   class   : vector of class labels corresponding to glomeruli (optional)
%
% OUTPUTS
%   distances : distance-to-kidney-edge at each glomerulus centroid
%   volumes   : volume of each glomerulus (in physical units)
%   cellcount : number of distinct cells per glomerulus (if inputs provided)
%   indices   : indices of glomeruli that were removed during filtering

% Voxel size
voxelSize = 4;  % assuming isotropic 4x4x4

% Measure glomerulus properties
props = regionprops3(volglom, 'Centroid', 'Volume');

% Number of glomeruli
numGlom = size(props, 1);

% Centroids returned as [x y z]; convert to integer voxel coordinates
centroids = round(cat(1, props.Centroid));   % [N x 3]

% Extract volumes (in voxels)
volumes = props.Volume;

% Convert centroid coordinates to linear indices
linIdx = sub2ind(size(volglom), ...
                 centroids(:,2), ... % row (y)
                 centroids(:,1), ... % column (x)
                 centroids(:,3));    % slice (z)

% Identify invalid indices (should be rare)
badIdx = isnan(linIdx);
indices = find(badIdx);
if ~isempty(indices)
    disp('Problem indices detected (NaN centroids)')
end

% Remove invalid entries
linIdx(badIdx) = [];
volumes(badIdx) = [];

% Distance to kidney edge at each centroid
distances = distK(linIdx);

% --- Optional: count distinct cells per glomerulus (VECTORIZED) ---
cellcount = [];

if exist('volcell', 'var') && ~isempty(volcell) && ...
   exist('D', 'var') && ~isempty(D) && ...
   exist('class', 'var') && ~isempty(class)

    % Ensure volcell is double for processing
    volcell = double(volcell);
    
    % Create mask where both glomerulus and cell labels are valid
    validMask = (volglom(:) > 0) & (volcell(:) > 0);
    
    % Extract glomerulus IDs and cell IDs at valid locations
    cellIDs = volglom(validMask);
    
    
    % Count distinct cells per glomerulus using accumarray
    rawCounts = accumarray(cellIDs, 1, [max(cellIDs), 1]);
    
    % Apply correction factor based on cell diameter
    correctionFactor = voxelSize / (voxelSize + mean(D(class)));
    cellcount = rawCounts * correctionFactor;
    
    % Remove bad indices from cellcount
    cellcount(badIdx) = [];
end

% --- Remove empty glomeruli ---
removeIdx = find(volumes == 0);

volumes(removeIdx)   = [];
distances(removeIdx) = [];

if ~isempty(cellcount)
    cellcount(removeIdx) = [];
end

indices = union(indices, removeIdx);

% --- Remove abnormally large glomeruli (volume threshold) ---
% Convert voxel count to physical volume
physicalVolumes = volumes * voxelSize^3;
largeIdx = find(physicalVolumes > 700000);

volumes(largeIdx)   = [];
distances(largeIdx) = [];

if ~isempty(cellcount)
    cellcount(largeIdx) = [];
end

indices = union(indices, largeIdx);

% Convert volumes to physical units
volumes = volumes * voxelSize^3;

end