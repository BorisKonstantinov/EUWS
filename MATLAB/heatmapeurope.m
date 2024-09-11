% figure;
% 
% % Geographical extents
% geoExtents = [-12, 32, 42, 72];
% 
% % Topographical map settings
% load topo;
% topoRefVec = [1, 90, 0];
% 
% latlim = [42 72];
% lonlim = [348 32];
% 
% worldmap(latlim, lonlim);
% geoshow(topo, topoRefVec, 'DisplayType', 'texturemap');
% colormap(topomap1);
% 
% title('Topographical Map');
% Assuming you have a matrix X with the wind data


list = generateList(energyHeatmapNAOPos);
list = list(list(:,3) ~= 0, :);
list = list(~isnan(list(:,3)), :);
list(:,3) = list(:,3)/10^9;

function result = generateList(X)
    % Check if the matrix X has the correct size
    if size(X) ~= [301, 471]
        error('Matrix X does not have the required size of 301x471.');
    end

    % Initialize the starting coordinates
    start_lat = 72.0;
    start_lon = -12.0;
    cell_size = 0.1;

    % Initialize the result array
    result = zeros(301 * 471, 3);  % Preallocate for speed

    % Fill the result array
    row_counter = 1;  % To keep track of rows in the result array

    for i = 1:301
        for j = 1:471
            result(row_counter, 1) = start_lat - (i-1) * cell_size;  % Latitude
            result(row_counter, 2) = start_lon + (j-1) * cell_size;  % Longitude
            result(row_counter, 3) = X(i, j);  % WindSpeed
            row_counter = row_counter + 1;
        end
    end
end

