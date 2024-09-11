figure;

% Geographical extents
geoExtents = [-12, 32, 42, 72];

% Create grid for your heatmap data
[longrid, latgrid] = meshgrid(linspace(geoExtents(1), geoExtents(2), size(energyHeatmap, 2)), ...
                             linspace(geoExtents(3), geoExtents(4), size(energyHeatmap, 1)));

% Plot the heatmap data
imagesc(longrid(1, :), latgrid(:, 1), energyHeatmap);
axis xy; % Ensure the data is not flipped vertically

% Custom colormap from white to deep red
nColors = 256; % Typically 256 for colormap
customColormap = [linspace(1, 0.5, nColors)', linspace(1, 0, nColors)', linspace(1, 0, nColors)'];
colormap(customColormap);
colorbar;

title('Heatmap');