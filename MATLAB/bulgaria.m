source_dir = 'C:/Users/bibob/Downloads';
target_dir = '2mtemp_Sofia_1940_2023.nc';
datafile = fullfile(source_dir, target_dir);


% Load latitude, longitude, time and temperature data
latitude = ncread(datafile, 'lat');
longitude = ncread(datafile, 'lon');
time = ncread(datafile, 'time');
temperature = ncread(datafile, 't2m');

% Convert the time from seconds since 1970 to datetime
time = datetime(time, 'ConvertFrom', 'posixtime');

% Convert the temperature from Kelvin to Celsius
temperature = temperature - 273.15;

% Select temperature data at the first longitude point
temperature_at_first_longitude = squeeze(temperature(1, 1, :));

% Extract indices of 04:00 and 16:00 data
indices_04 = (hour(time) == 4);
indices_16 = (hour(time) == 16);

% Extract 04:00 and 16:00 data
time_04 = time(indices_04);
temperature_04 = temperature_at_first_longitude(indices_04);

time_16 = time(indices_16);
temperature_16 = temperature_at_first_longitude(indices_16);

% Set the size of the scatter plot markers
markerSize = 10;

% Fit a trendline and calculate percentiles for 04:00 data
p_04 = polyfit(datenum(time_04), temperature_04, 1);
temp_04_25 = prctile(temperature_04, 2.5);
temp_04_75 = prctile(temperature_04, 97.5);

% Create scatter plot for 04:00 data
figure;
scatter(time_04, temperature_04, markerSize, 'b');
hold on;
plot(time_04, polyval(p_04, datenum(time_04)), '-g', 'LineWidth', 2); % trendline
line(get(gca, 'XLim'), [temp_04_25 temp_04_25], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % 25th percentile
line(get(gca, 'XLim'), [temp_04_75 temp_04_75], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % 75th percentile
hold off;
xlabel('Time');
ylabel('Temperature (°C)');
title('Temperature at Sofia from 1940 to 2023 at 04:00');
grid on;
set(gca, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5, 'YGrid', 'on', 'YMinorGrid', 'off');
set(gca, 'YTick', floor(min(temperature_04)):1:ceil(max(temperature_04)));

% % Draw thicker black lines every 5°C for 04:00 data
% for y = floor(min(temperature_04)):5:ceil(max(temperature_04))
%     line(get(gca, 'XLim'), [y y], 'Color', 'k', 'LineWidth', 1.5);
% end

% Fit a trendline and calculate percentiles for 16:00 data
p_16 = polyfit(datenum(time_16), temperature_16, 1);
temp_16_25 = prctile(temperature_16, 2.5);
temp_16_75 = prctile(temperature_16, 97.5);

% Create scatter plot for 16:00 data
figure;
scatter(time_16, temperature_16, markerSize, 'r');
hold on;
plot(time_16, polyval(p_16, datenum(time_16)), '-g', 'LineWidth', 2); % trendline
line(get(gca, 'XLim'), [temp_16_25 temp_16_25], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % 25th percentile
line(get(gca, 'XLim'), [temp_16_75 temp_16_75], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % 75th percentile
hold off;
xlabel('Time');
ylabel('Temperature (°C)');
title('Temperature at Sofia from 1940 to 2023 at 16:00');
grid on;
set(gca, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5, 'YGrid', 'on', 'YMinorGrid', 'off');
set(gca, 'YTick', floor(min(temperature_16)):1:ceil(max(temperature_16)));

% % Draw thicker black lines every 5°C for 16:00 data
% for y = floor(min(temperature_16)):5:ceil(max(temperature_16))
%     line(get(gca, 'XLim'), [y y], 'Color', 'k', 'LineWidth', 1.5);
% end