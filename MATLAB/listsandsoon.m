% listPminusNeg = subtractLists(listNAOPositive, listNAONegative);
% listPminusNeu = subtractLists(listNAOPositive, listNAONegative);
% listNeuminusNeg = subtractLists(listNAOPositive, listNAONegative);
% 
% function result = subtractLists(list1, list2)
%     % Combining both lists to ensure we capture all unique coordinates
%     combined_coords = [list1(:, 1:2); list2(:, 1:2)];
%     [unique_coords, ~, ic] = unique(combined_coords, 'rows', 'stable');
% 
%     % Initialize result array
%     result = [unique_coords, zeros(size(unique_coords, 1), 1)];
% 
%     for i = 1:size(unique_coords, 1)
%         % Extracting wind speeds for the current coordinates from both lists
%         ind1 = ismember(list1(:, 1:2), unique_coords(i, :), 'rows');
%         ind2 = ismember(list2(:, 1:2), unique_coords(i, :), 'rows');
% 
%         ws1 = 0;
%         ws2 = 0;
% 
%         if any(ind1)
%             ws1 = list1(ind1, 3);
%         end
% 
%         if any(ind2)
%             ws2 = list2(ind2, 3);
%         end
% 
%         % Subtracting wind speeds
%         result(i, 3) = ws1 - ws2;
%     end
% end




% Create a new figure
figure;

% Generate a colormap that smoothly transitions between the two colors
numColors = 256; % Number of steps in the gradient
colorStart = [255, 255, 102] / 255;
colorEnd = [255, 0, 0] / 255;

R = linspace(colorStart(1), colorEnd(1), numColors)';
G = linspace(colorStart(2), colorEnd(2), numColors)';
B = linspace(colorStart(3), colorEnd(3), numColors)';

customColormap = [R, G, B];
colormap(customColormap);

% Create a dummy plot (this is needed to then display a colorbar)
imagesc([0 1; 0 1]);

% Display the colorbar
c = colorbar;
c.Label.String = 'Giga Joules (GJ)';
caxis([0 10000000]); % Set the limits for the colorbar

% Hide the dummy plot
axis off;

% Save the colorbar as an image
saveas(gcf, 'colorbar_gradient.png');

