% Figure 1: Plot daily index and monthly index for target year
figure(1);
clf;
% Specify target year
targetYear = 2020;
rows = year(naoDpM.date) == targetYear;
% 
bar(naoDpM.date(rows), naoDpM.index_d(rows), 'FaceColor', [0.7 0.7 0.7], 'DisplayName', 'Daily Index');
hold on;
plot(naoDpM.date(rows), naoDpM.index_m(rows), 'DisplayName', 'Monthly Index');
yline(0, 'DisplayName', 'y=0');
hold off;
% Labeling
xticks(datetime(targetYear,1:12,1));
grid on
title(sprintf('Index values for year: %d', targetYear));
legend('show');
xlabel('Date');
ylabel('Index Value');

% Figure 2: Scatter plot of daily index vs monthly index
figure(2);
clf;
scatter(naoDpM.index_d, naoDpM.index_m);
    hold on      
yline(0, 'DisplayName', 'y=0');
hold on
xline(0, 'DisplayName', 'x=0');
% Labeling
grid on;
legend('show');
title('Daily Index vs Monthly Index');
xlabel('Daily Index');
ylabel('Monthly Index');
hold off;

% Figure 3: Scatter plot of average daily index vs monthly index with best fit line
figure(3);
clf;
scatter(naoDpM.index_d_avg, naoDpM.index_m);
hold on;
% Compute and plot best fit line
p = polyfit(naoDpM.index_d_avg, naoDpM.index_m, 1); % Fit a first degree polynomial
yfit = polyval(p, naoDpM.index_d_avg); % Get y-values of the best fit line
plot(naoDpM.index_d_avg, yfit, 'r-', 'DisplayName', 'Best fit line'); % Plot the best fit line
% Create equation textbox
equationText = sprintf('Best fit line: y = %.3f*x + %.3f', p(1), p(2));
dim = [.6 .2 .3 .3]; % Adjust these values to move the textbox position
annotation('textbox',dim,'String',equationText,'FitBoxToText','on');
% Labeling
title('Average Daily Index vs Monthly Index');
xlabel('Average Daily Index');
ylabel('Monthly Index');
grid on;
legend('show');
hold off;

% Daily max wind speed in 30:70N, -80:20E
figure(4);
clf;

% Specify target year
targetYear = 2020;
rows = year(naoDpM.date) == targetYear;

% Find the maximum daily NAO value for the target year
maxNAOValue = max(naoDpM.index_d(rows));

% Find the maximum wind speed in the table
maxWindSpeed = max(dailyMaxWindSpeed.max_wind_speed);

% Plot daily index and monthly index
bar(naoDpM.date(rows), naoDpM.index_d(rows), 'FaceColor', [0.7 0.7 0.7], 'DisplayName', 'Daily Index');
hold on;
plot(naoDpM.date(rows), naoDpM.index_m(rows), 'DisplayName', 'Monthly Index');

% Plot daily maximum wind speed
% Extract the data for the target year
windRows = year(dailyMaxWindSpeed.date) == targetYear;

% Scale the wind speed by its maximum value, then by the maximum daily NAO value
scaledWindSpeed = (dailyMaxWindSpeed.max_wind_speed(windRows) / maxWindSpeed) * maxNAOValue;

% Plot the daily maximum wind speed
plot(dailyMaxWindSpeed.date(windRows), scaledWindSpeed, 'o', 'DisplayName', 'Daily Max Wind Speed');

% Add reference line at y=0
yline(0, 'DisplayName', 'y=0');

hold off;

% Labeling
xticks(datetime(targetYear,1:12,1));
grid on
title(sprintf('Index values and wind speed for year: %d', targetYear));
legend('show');
xlabel('Date');
ylabel('Index Value / Wind Speed');
