% Figure1(naoDpM, 2020);
% plot_nao_yearly_avg(naoDpM);
% Figure2(naoDpM);
% Figure3(naoDpM);
% Figure4(naoDpM, 2020, EUWS);
% Figure5(SE);
% Figure6(SESSI);
% Figure6(SE_BEL);
% Figure6(SE_ESP);
% Figure7(SESSI, naoDpM)

% plotNAOIndexChangeDuringStorms(SESSI, naoDpM);
% plot_cum_prob(SESSI, naoDpM);

% plot_nao_windstorm(naoDpM, SESSI);

% plotSESSI(SESSI);

% Figure6b(SESSI, 2000, 2020);

% plot_windstorm_energy_fraction(SESSI, SE_AUT, SE_BEL, SE_CHE, SE_CZE, SE_DEU, SE_DNK, SE_ESP, SE_EST, SE_FIN, SE_FRA, SE_GBR, SE_HUN, SE_IRL, SE_LTU, SE_LVA, SE_NLD, SE_NOR, SE_POL, SE_SVK, SE_SWE);
% plot_average_windstorm_energy_fraction(SESSI, SE_AUT, SE_BEL, SE_CHE, SE_CZE, SE_DEU, SE_DNK, SE_ESP, SE_EST, SE_FIN, SE_FRA, SE_GBR, SE_HUN, SE_IRL, SE_LTU, SE_LVA, SE_NLD, SE_NOR, SE_POL, SE_SVK, SE_SWE);
% plot_total_windstorm_energy(SESSI, SE_AUT, SE_BEL, SE_CHE, SE_CZE, SE_DEU, SE_DNK, SE_ESP, SE_EST, SE_FIN, SE_FRA, SE_GBR, SE_HUN, SE_IRL, SE_LTU, SE_LVA, SE_NLD, SE_NOR, SE_POL, SE_SVK, SE_SWE);
% plot_windstorm_event_count(SESSI, SE_AUT, SE_BEL, SE_CHE, SE_CZE, SE_DEU, SE_DNK, SE_ESP, SE_EST, SE_FIN, SE_FRA, SE_GBR, SE_HUN, SE_IRL, SE_LTU, SE_LVA, SE_NLD, SE_NOR, SE_POL, SE_SVK, SE_SWE);


%Figure6(SE_ESP);
%Figure6(SE_FRA);
%Figure6(SE_BEL);
%Figure6(SE_NLD);
%Figure6(SE_GBR);
%Figure6(SE_DEU);
%Figure6(SE_DNK);
%Figure6(SE_NOR);





% Daily + Monthly index over time
function Figure1(naoDpM, targetYear)
figure(1);
clf;
% Setup
rows = year(naoDpM.date) == targetYear;
% Plotting
plot(naoDpM.date(rows), naoDpM.index_d(rows), '--.k', 'DisplayName', 'Daily Index');
hold on;
plot(naoDpM.date(rows), naoDpM.index_m(rows), '-b', 'DisplayName', 'Monthly Index', 'LineWidth', 1.2);
% Readability
title(sprintf('Index values for year: %d', targetYear));
xlabel('Date');
ylabel('Index Value');
grid on;
legend('show');
xticks(datetime(targetYear,1:12,1));
yline(0, 'DisplayName', 'y=0');

hold off;
end

% Scatter plot of Daily vs Monthly index
function Figure2(naoDpM)
figure(2);
clf;
% Setup
% Plotting
scatter(naoDpM.index_d, naoDpM.index_m);
hold on
yline(0, 'DisplayName', 'y=0');
hold on
xline(0, 'DisplayName', 'x=0');
% Readability
title('Daily Index vs Monthly Index');
xlabel('Daily Index');
ylabel('Monthly Index');
grid on;
legend('show');

hold off;
end

% Scatter plot of Avg Daily vs Monthly index + line best fit
function Figure3(naoDpM)
figure(3);
clf;
% Setup
p = polyfit(naoDpM.index_d_avg, naoDpM.index_m, 1); % Fit a first degree polynomial
yfit = polyval(p, naoDpM.index_d_avg); % Get y-values of the best fit line
equationText = sprintf('Best fit line: y = %.3f*x + %.3f', p(1), p(2));
% Plotting
scatter(naoDpM.index_d_avg, naoDpM.index_m);
hold on;
plot(naoDpM.index_d_avg, yfit, 'r-', 'DisplayName', 'Best fit line'); % Plot the best fit line
% Readability
title('Average Daily Index vs Monthly Index');
xlabel('Average Daily Index');
ylabel('Monthly Index');
grid on;
legend('show');
annotation('textbox', [.6 .2 .3 .3], 'String', equationText, 'FitBoxToText','on');

hold off;
end

% Daily max wind speed in 30:70N, -80:20E (WTF the point of this bich?)
function Figure4(naoDpM, targetYear, EUWS)
figure(4);
clf;
% Setup
rows = year(naoDpM.date) == targetYear;
maxNAOValue = max(naoDpM.index_d(rows));
dailyMaxWindSpeed = varfun(@max, EUWS, 'GroupingVariables', 'date', 'InputVariables', 'wind_speed_10m');
maxWindSpeed = max(dailyMaxWindSpeed.max_wind_speed_10m);
windRows = year(dailyMaxWindSpeed.date) == targetYear;
scaledWindSpeed = (dailyMaxWindSpeed.max_wind_speed_10m(windRows) / maxWindSpeed) * maxNAOValue;
% Plotting
bar(naoDpM.date(rows), naoDpM.index_d(rows), 'FaceColor', [0.7 0.7 0.7], 'DisplayName', 'Daily Index');
hold on;
plot(naoDpM.date(rows), naoDpM.index_m(rows), 'DisplayName', 'Monthly Index');
hold on;
plot(dailyMaxWindSpeed.date(windRows), scaledWindSpeed, 'o', 'DisplayName', 'Daily Max Wind Speed');
% Readability
title(sprintf('Index values and wind speed for year: %d', targetYear));
xlabel('Date');
ylabel('Index Value / Wind Speed');
grid on
legend('show');
xticks(datetime(targetYear,1:12,1));
yline(0, 'DisplayName', 'y=0');

hold off;
end

% Plots Storm Energy for a period Bar style with Width = Storm Length
function Figure5(SE)
figure(5);
clf;
% Setup
SE = SE(1:20, :);
stormDuration = SE.endTime - SE.startTime;
stormDurationDays = days(stormDuration);
% Plotting
for i = 1:height(SE)
    % Compute the x, y, width, and height of the rectangle to be drawn
    x = SE.startTime(i);
    y = 0;
    width = stormDuration(i);
    heightSE = SE.stormEnergy(i);
    % Plot a rectangle (patch) for this row
    patch([x x x+width x+width], [y y+heightSE y+heightSE y], 'b');
    hold on;
end
% Readability
title('Storm Energy over Time');
xlabel('Time');
ylabel('Storm Energy');
grid on
legend('show');
end

% Plots Exceedance Probability curve for NAO(+, 0, -) for 1950->2020
function Figure6(SESSI)
    figure();
    % Setup
    positiveNAO = SESSI.stormEnergy(SESSI.NAOIndex > 0.5 & SESSI.stormEnergy ~= 0);
    neutralNAO = SESSI.stormEnergy(SESSI.NAOIndex <= 0.5 & SESSI.NAOIndex >= -0.5 & SESSI.stormEnergy ~= 0);
    negativeNAO = SESSI.stormEnergy(SESSI.NAOIndex < -0.5 & SESSI.stormEnergy ~= 0);

    % Remove NaN values
    positiveNAO = positiveNAO(~isnan(positiveNAO));
    neutralNAO = neutralNAO(~isnan(neutralNAO));
    negativeNAO = negativeNAO(~isnan(negativeNAO));

    % Calculate EP
    N_positive = length(positiveNAO);
    N_neutral = length(neutralNAO);
    N_negative = length(negativeNAO);

    rank_positive = (1:N_positive)';
    rank_neutral = (1:N_neutral)';
    rank_negative = (1:N_negative)';

    P_positive = rank_positive / (N_positive + 1); % Weibull plotting position
    P_neutral = rank_neutral / (N_neutral + 1);    % Weibull plotting position
    P_negative = rank_negative / (N_negative + 1); % Weibull plotting position

    % Sort data in descending order
    sorted_positive = sort(positiveNAO, 'descend');
    sorted_neutral = sort(neutralNAO, 'descend');
    sorted_negative = sort(negativeNAO, 'descend');

    % Plotting
    semilogx(sorted_positive, P_positive);
    hold on;
    semilogx(sorted_neutral, P_neutral);
    semilogx(sorted_negative, P_negative);
    hold off;

    % Readability
    title('Exceedance Probability Curve for Positive, Neutral, and Negative NAO Indices');
    xlabel('Storm Energy [unit]');
    ylabel('Exceedance Probability');

    % Compute stats and add to the legend
    legendStr = {
        ['Positive NAO,   N = ' num2str(N_positive) ', Mean: ' sprintf('%.1e', mean(positiveNAO)) '[J], Median: ' sprintf('%.1e', median(positiveNAO)) '[J]'],
        ['Neutral NAO,     N = ' num2str(N_neutral) ', Mean: ' sprintf('%.1e', mean(neutralNAO)) '[J], Median: ' sprintf('%.1e', median(neutralNAO)) '[J]'],
        ['Negative NAO, N = ' num2str(N_negative) ', Mean: ' sprintf('%.1e', mean(negativeNAO)) '[J], Median: ' sprintf('%.1e', median(negativeNAO)) '[J]']
        };
    legend(legendStr, 'Location', 'best');
end

% Create Bar plot of the storm count during NAO(+, 0, -) for 1950->2020
function Figure7(SESSI, naoDpM)
figure(7);
clf
% Setup
edges = -3.5:0.5:3.5; % Edit as necessary
stormCountMonthly = zeros(length(edges)-1, 1);
stormCountDaily = zeros(length(edges)-1, 1);
for i = 1:length(edges)-1
    % Define the current bin.
    bin = [edges(i) edges(i+1)];

    % Find the storms that fall within the current bin, for the monthly index.
    idx = (naoDpM.index_m >= bin(1)) & (naoDpM.index_m < bin(2));
    relevantNAOMonths = naoDpM.date(idx);
    stormCountMonthly(i) = sum(ismember(floor(datenum(SESSI.startTime)), datenum(relevantNAOMonths)));

    % Repeat for the daily index.
    idx = (naoDpM.index_d >= bin(1)) & (naoDpM.index_d < bin(2));
    relevantNAODays = naoDpM.date(idx);
    stormCountDaily(i) = sum(ismember(floor(datenum(SESSI.startTime)), datenum(relevantNAODays)));
end
yMax = max(max(stormCountMonthly), max(stormCountDaily));
xMin = edges(1);
xMax = edges(end);
% Plotting Subplot 1
subplot(2, 1, 1);
bar(edges(1:end-1), stormCountMonthly, 'histc');
xlim([xMin xMax]);
ylim([0 yMax]);
% Readability Subplot 1
title('Storm Count vs. Monthly NAO Index');
xlabel('NAO Index');
ylabel('Number of Storms');
% Plotting Subplot 2
subplot(2, 1, 2);
bar(edges(1:end-1), stormCountDaily, 'histc');
xlim([xMin xMax]);
ylim([0 yMax]);
% Readability Subplot 2
title('Storm Count vs. Daily NAO Index');
xlabel('NAO Index');
ylabel('Number of Storms');
end

% 
function plotNAOIndexChangeDuringStorms(SESSI, naoDpM)
    % Initialize arrays to store the NAO index at the start and end of each storm.
    startNAOIndex = nan(height(SESSI), 1);
    endNAOIndex = nan(height(SESSI), 1);

    % Loop over the storms.
    for i = 1:height(SESSI)
        % Get the start and end times of the storm.
        startTime = SESSI.startTime(i);
        endTime = SESSI.endTime(i);

        % Find the NAO indices corresponding to the start and end times.
        startIndex = find(naoDpM.date == dateshift(startTime,'start','day'));
        endIndex = find(naoDpM.date == dateshift(endTime,'start','day'));
        
        if ~isempty(startIndex)
            startNAOIndex(i) = naoDpM.index_d(startIndex(1));
        end
        if ~isempty(endIndex)
            endNAOIndex(i) = naoDpM.index_d(endIndex(1));
        end
    end

    % Remove any entries where either the start or end NAO index was not found.
    validData = ~isnan(startNAOIndex) & ~isnan(endNAOIndex);
    startNAOIndex = startNAOIndex(validData);
    endNAOIndex = endNAOIndex(validData);

    % Plot the start NAO index against the end NAO index.
    figure;
    scatter(startNAOIndex, endNAOIndex);
    hold on;
    
    % Fit a straight line (1st degree polynomial) to the data
    p = polyfit(startNAOIndex, endNAOIndex, 1);
    
    % Generate the x and y values for the fitted line
    x_fit = linspace(min(startNAOIndex), max(startNAOIndex), 100);
    y_fit = polyval(p, x_fit);
    
    % Add the fitted line to the plot
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    
    % Display the equation of the best fit line
    equationText = sprintf('y = %.2fx + %.2f', p(1), p(2));
    text(min(startNAOIndex), max(endNAOIndex), equationText, 'FontSize', 12, 'Color', 'r');
    
    hold off;
    title('NAO Index at Start vs End of Storms');
    xlabel('NAO Index at Start');
    ylabel('NAO Index at End');
    legend('Data', 'Best fit line', 'Location', 'best');
end

%
function plot_cum_prob(SESSI, naoDpM)

% Convert the datetime to date
SESSI.startTime = dateshift(SESSI.startTime, 'start', 'day');
SESSI.endTime = dateshift(SESSI.endTime, 'start', 'day');

% Remove rows with NaN values in 'index_d'
naoDpM = rmmissing(naoDpM,'DataVariables',{'index_d'});

% Classify NAO states based on 'index_d'
naoDpM.naoState(naoDpM.index_d > 0.5) = {'positive'};
naoDpM.naoState(naoDpM.index_d <= 0.5 & naoDpM.index_d >= -0.5) = {'neutral'};
naoDpM.naoState(naoDpM.index_d < -0.5) = {'negative'};

% Create a timetable for naoDpM
naoDpM_tt = table2timetable(naoDpM);

% Calculate total number of storms
totalStorms = height(SESSI);

% Count the number of storms each month
stormsPerMonth = histcounts(month(SESSI.startTime), 1:13);

% Calculate storm probabilities for each month
stormProb = stormsPerMonth / totalStorms;

% Initialize storm counts matrix (rows: NAO states, columns: months)
stormCounts = zeros(3,12);

% Loop over the SESSI table
for i = 1:height(SESSI)
    stormStart = SESSI.startTime(i);
    
    % Select the relevant NAO state
    relevantNaoState = naoDpM_tt(naoDpM_tt.date == stormStart, :).naoState;
    
    % Update the storm counts based on the NAO state and the month of the storm
    if ~isempty(relevantNaoState)
        stormMonth = month(stormStart);

        if strcmp(relevantNaoState, 'positive')
            stormCounts(1, stormMonth) = stormCounts(1, stormMonth) + 1;
        elseif strcmp(relevantNaoState, 'neutral')
            stormCounts(2, stormMonth) = stormCounts(2, stormMonth) + 1;
        else
            stormCounts(3, stormMonth) = stormCounts(3, stormMonth) + 1;
        end
    end
end

% Calculate the probabilities for each NAO state given the month and normalize them
stormProb_NAO = stormCounts / totalStorms;

% Create a bar chart for overall storm probabilities
figure;
bar(1:12, stormProb);
xlabel('Month');
ylabel('Probability');
title('Probability of Storm Events');

% Create a bar chart for normalized probabilities given NAO state
figure;
bar(1:12, stormProb_NAO', 'stacked');
xlabel('Month');
ylabel('Normalized Probability');
legend('NAO Positive', 'NAO Neutral', 'NAO Negative');
title('Normalized Probability of Storm Events Given NAO State');

end

%
function plot_nao_windstorm(naoDpM, SESSI)

    % Define NAO index bins
    nao_bins = -3.5:0.5:3.5;
    nao_bin_midpoints = nao_bins(1:end-1) + 0.5 * diff(nao_bins);

    % Extract NAO data and dates
    nao_daily = naoDpM.index_d;
    nao_monthly = naoDpM.index_m;
    nao_dates = datevec(naoDpM.date);

    % Extract event dates
    event_dates = datevec(SESSI.startTime);

    % Create figure
    figure;

    % Subplot 1: Daily NAO histogram
    subplot(2,2,1);
    histogram(nao_daily, nao_bins);
    title('Daily NAO distribution');
    xlabel('NAO Index');
    ylabel('Day count');

    % Subplot 2: Event ratios for daily NAO
    subplot(2,2,2);
    event_ratios_daily = calculate_event_ratios(nao_dates, nao_daily, event_dates, nao_bin_midpoints);
    bar(nao_bin_midpoints, event_ratios_daily);
    title('Event ratios (Daily NAO)');
    xlabel('NAO Index');
    ylabel('Event ratio');

    % Subplot 3: Monthly NAO histogram
    subplot(2,2,3);
    histogram(nao_monthly, nao_bins);
    title('Monthly NAO distribution');
    xlabel('NAO Index');
    ylabel('Month count');

    % Subplot 4: Event ratios for monthly NAO
    subplot(2,2,4);
    event_ratios_monthly = calculate_event_ratios(nao_dates, nao_monthly, event_dates, nao_bin_midpoints);
    bar(nao_bin_midpoints, event_ratios_monthly);
    title('Event ratios (Monthly NAO)');
    xlabel('NAO Index');
    ylabel('Event ratio');
end

%
function event_ratios = calculate_event_ratios(nao_dates, nao_indices, event_dates, nao_bin_midpoints)

    % Preallocate event ratios
    event_ratios = zeros(size(nao_bin_midpoints));

    % For each NAO bin, calculate the ratio of events to days/months
    for i = 1:numel(nao_bin_midpoints)

        % Find days/months with NAO in this bin
        nao_in_bin = abs(nao_indices - nao_bin_midpoints(i)) < 0.25;

        % Find events on these days/months
        events_in_bin = ismember(event_dates(:,1:3), nao_dates(nao_in_bin,1:3), 'rows');

        % Calculate event ratio
        event_ratios(i) = sum(events_in_bin) / sum(nao_in_bin);
    end
end

%
function plot_nao_yearly_avg(data)
    % Assume data is a table with columns: index_d, index_m, date

    % Convert dates to datetime format if not already
    if ischar(data.date(1)) || isstring(data.date(1))
        data.date = datetime(data.date, 'InputFormat', 'dd-MMM-yyyy');
    end
    
    % Extract years
    years = year(data.date);
    
    % Calculate yearly average using the monthly index
    [unique_years, ~, idy] = unique(years);
    yearly_avg = accumarray(idy, data.index_m, [], @mean);
    
    % Generate colors based on values
    peak_intensity = 2.0;
    colors = zeros(length(yearly_avg), 3);
    for i = 1:length(yearly_avg)
        if yearly_avg(i) >= 0
            intensity = min(yearly_avg(i) / peak_intensity, 1);
            colors(i, :) = [1, 1-intensity, 1-intensity];
        else
            intensity = min(-yearly_avg(i) / peak_intensity, 1);
            colors(i, :) = [1-intensity, 1-intensity, 1];
        end
    end
    
    % Plot
    figure;
    b = bar(unique_years, yearly_avg, 'FaceColor', 'flat', 'EdgeColor', 'k');
    for i = 1:length(yearly_avg)
        b.CData(i, :) = colors(i, :);
    end
    title('Average Yearly NAO Index (Using Monthly Data)');
    xlabel('Year');
    ylabel('Average Yearly NAO Index');
    grid on;

    % Add a colorbar
    c = colorbar;
    caxis([-peak_intensity, peak_intensity]);  % Set colorbar limits
    ylabel(c, 'NAO Index Value');
    c.Ticks = [-peak_intensity, 0, peak_intensity];
    c.TickLabels = {sprintf('%.1f', -peak_intensity), '0', sprintf('%.1f', peak_intensity)};
end

%
function plotSESSI(SESSI)
    % Extract month and categorize NAO index
    months = month(datetime(SESSI.startTime, 'InputFormat', 'yyyy-MM-dd HH:mm:ss'));
    bins = [-inf, -0.5, 0.5, inf];
    labels = {'Negative', 'Neutral', 'Positive'};
    NAO_category = discretize(SESSI.NAOIndex, bins, 'Categorical', labels);

    % Initialize data for plotting and counting
    barData = zeros(12, 3);
    counts = zeros(12, 3);  % Store counts of events

    for i = 1:12
        for j = 1:3
            idx = (months == i) & (NAO_category == labels{j});
            counts(i, j) = sum(idx);
            if any(idx)
                barData(i, j) = nanmean(SESSI.stormEnergy(idx));  % Use nanmean to ignore NaN values
            end
        end
    end

    % Plot
    figure;
    b = bar(1:12, barData, 'grouped');
    title('Average Storm Energy by Month for Different NAO Categories', 'FontSize', 15, 'FontName', 'Courier New', 'FontWeight', 'bold');
    xlabel('Month', 'FontSize', 15, 'FontName', 'Courier New', 'FontWeight', 'bold');
    ylabel('Average Storm Energy', 'FontSize', 15, 'FontName', 'Courier New', 'FontWeight', 'bold');
    legend(labels, 'Location', 'northwest', 'FontSize', 15, 'FontName', 'Courier New', 'FontWeight', 'bold');

    % Adjust x-ticks for clarity
    xticks(1:12);
    xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'});
    set(gca, 'FontSize', 15, 'FontName', 'Courier New', 'FontWeight', 'bold');

    % Annotate bars with counts
    for k = 1:12
        for m = 1:3
            height = barData(k, m);
            if height > 0
                count = counts(k, m);
                barCenter = b(m).XEndPoints(k); % Get the center of the current bar
                text(barCenter, 0, sprintf('N = %d', count), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k', 'FontSize', 15, 'FontName', 'Courier New', 'Rotation', 90, 'FontWeight', 'bold');
            end
        end
    end
end

%
function Figure6b(SESSI, startYear, endYear)
    figure();
    
    % Filter by year
    yearFilter = year(SESSI.startTime) >= startYear & year(SESSI.startTime) <= endYear;

    % Setup
    positiveNAO = SESSI.stormEnergy(yearFilter & SESSI.NAOIndex > 0.5 & SESSI.stormEnergy ~= 0);
    neutralNAO = SESSI.stormEnergy(yearFilter & SESSI.NAOIndex <= 0.5 & SESSI.NAOIndex >= -0.5 & SESSI.stormEnergy ~= 0);
    negativeNAO = SESSI.stormEnergy(yearFilter & SESSI.NAOIndex < -0.5 & SESSI.stormEnergy ~= 0);

    % Remove NaN values
    positiveNAO = positiveNAO(~isnan(positiveNAO));
    neutralNAO = neutralNAO(~isnan(neutralNAO));
    negativeNAO = negativeNAO(~isnan(negativeNAO));

    % Calculate EP
    N_positive = length(positiveNAO);
    N_neutral = length(neutralNAO);
    N_negative = length(negativeNAO);

    rank_positive = (1:N_positive)';
    rank_neutral = (1:N_neutral)';
    rank_negative = (1:N_negative)';

    P_positive = rank_positive / (N_positive + 1); % Weibull plotting position
    P_neutral = rank_neutral / (N_neutral + 1);    % Weibull plotting position
    P_negative = rank_negative / (N_negative + 1); % Weibull plotting position

    % Sort data in descending order
    sorted_positive = sort(positiveNAO, 'descend');
    sorted_neutral = sort(neutralNAO, 'descend');
    sorted_negative = sort(negativeNAO, 'descend');

    % Plotting
    semilogx(sorted_positive, P_positive);
    hold on;
    semilogx(sorted_neutral, P_neutral);
    semilogx(sorted_negative, P_negative);
    hold off;

    % Readability
    title(['Exceedance Probability Curve for Positive, Neutral, and Negative NAO Indices (' num2str(startYear) '-' num2str(endYear) ')']);
    xlabel('Storm Energy [unit]');
    ylabel('Exceedance Probability');

    % Compute stats and add to the legend
    legendStr = {
        ['Positive NAO,   N = ' num2str(N_positive) ', Mean: ' sprintf('%.1e', mean(positiveNAO)) '[J], Median: ' sprintf('%.1e', median(positiveNAO)) '[J]'],
        ['Neutral NAO,     N = ' num2str(N_neutral) ', Mean: ' sprintf('%.1e', mean(neutralNAO)) '[J], Median: ' sprintf('%.1e', median(neutralNAO)) '[J]'],
        ['Negative NAO, N = ' num2str(N_negative) ', Mean: ' sprintf('%.1e', mean(negativeNAO)) '[J], Median: ' sprintf('%.1e', median(negativeNAO)) '[J]']
        };
    legend(legendStr, 'Location', 'best');
end

%winter + summer
function plot_windstorm_energy_fraction(SESSI, varargin)
    % This function takes in the main SESSI table and country-specific tables 
    % as varargin and plots the fraction of windstorm energy for each country.
    
    % Filter out windstorm events that have their start time during December, January, or February
    winter_months = [6, 7, 8];
    SESSI.startTime = datetime(SESSI.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    SESSI_winter = SESSI(ismember(month(SESSI.startTime), winter_months),:);
    
    % Calculate the total energy of these windstorms
    total_energy_winter = sum(SESSI_winter.stormEnergy);
    
    % For each country, calculate the fraction of energy experienced by that country
    num_countries = length(varargin);
    countries = cell(1, num_countries);
    fractions = zeros(1, num_countries);
    for i = 1:num_countries
        country_table = varargin{i};
        country_table.startTime = datetime(country_table.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        country_winter = country_table(ismember(month(country_table.startTime), winter_months),:);
        fractions(i) = sum(country_winter.stormEnergy) / total_energy_winter;
        countries{i} = extractAfter(inputname(i+1), "SE_");  % Extract country code from variable name
    end
    
    % Plot
    figure;
    bar(fractions);
    ax = gca;
    ax.XTick = 1:num_countries; % Set XTick to position of each bar
    ax.XTickLabel = countries;
    set(gca, 'XTickLabelRotation', 45); % Rotate x-axis labels by 45 degrees for better visibility
    xlabel('Country');
    ylabel('Fraction of Windstorm Energy');
    title('Fraction of Boreal Winter Windstorm Energy Experienced by Country');
    ylim([0 1]);
    grid on;
end

% dont need
function plot_average_windstorm_energy_fraction(SESSI, varargin)
    % This function takes in the main SESSI table and country-specific tables 
    % as varargin and plots the average fraction of windstorm energy for each country.
    
    % Filter out windstorm events that have their start time during December, January, or February
    winter_months = [6, 7, 8];
    SESSI.startTime = datetime(SESSI.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    SESSI_winter = SESSI(ismember(month(SESSI.startTime), winter_months),:);
    
    % For each country, calculate the fraction of energy experienced by that country for each storm
    num_countries = length(varargin);
    num_storms = height(SESSI_winter);
    countries = cell(1, num_countries);
    fractions = zeros(num_countries, num_storms);
    for i = 1:num_countries
        country_table = varargin{i};
        country_table.startTime = datetime(country_table.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        for j = 1:num_storms
            storm_start_time = SESSI_winter.startTime(j);
            storm_energy = SESSI_winter.stormEnergy(j);
            
            % Find the corresponding entry in the country table based on startTime
            country_entry = country_table(country_table.startTime == storm_start_time, :);
            if isempty(country_entry)
                country_storm_energy = 0;
            else
                country_storm_energy = country_entry.stormEnergy;
            end
            
            fractions(i, j) = country_storm_energy / storm_energy;
        end
        countries{i} = extractAfter(inputname(i+1), "SE_");  % Extract country code from variable name
    end
    
    % Calculate the average fraction for each country
    average_fractions = nanmean(fractions, 2);  % Use nanmean to ignore NaN values
    
    % Plot
    figure;
    bar(average_fractions');
    ax = gca;
    ax.XTick = 1:num_countries;
    ax.XTickLabel = countries;
    set(gca, 'XTickLabelRotation', 45); % Rotate x-axis labels by 45 degrees for better visibility
    xlabel('Country');
    ylabel('Average Fraction of Windstorm Energy');
    title('Average Fraction of Boreal Winter Windstorm Energy Experienced by Country');
    ylim([0 1]);
    grid on;
end

%winter
function plot_total_windstorm_energy(SESSI, varargin)
    % This function takes in the main SESSI table and country-specific tables 
    % as varargin and plots the total windstorm energy for each country during boreal winter.
    
    % Filter out windstorm events that have their start time during December, January, or February
    winter_months = [6, 7, 8];
    SESSI.startTime = datetime(SESSI.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    SESSI_winter = SESSI(ismember(month(SESSI.startTime), winter_months),:);
    
    % For each country, sum the energy experienced by that country for each storm
    num_countries = length(varargin);
    countries = cell(1, num_countries);
    total_energies = zeros(1, num_countries);
    for i = 1:num_countries
        country_table = varargin{i};
        country_table.startTime = datetime(country_table.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        country_winter = country_table(ismember(month(country_table.startTime), winter_months),:);
        total_energies(i) = sum(country_winter.stormEnergy);
        countries{i} = extractAfter(inputname(i+1), "SE_");  % Extract country code from variable name
    end
    
    % Plot
    figure;
    bar(total_energies);
    ax = gca;
    ax.XTick = 1:num_countries;
    ax.XTickLabel = countries;
    set(gca, 'XTickLabelRotation', 45); % Rotate x-axis labels by 45 degrees for better visibility
    xlabel('Country');
    ylabel('Total Windstorm Energy');
    title('Total Boreal Winter Windstorm Energy Experienced by Country');
    grid on;
end

%winter + summer
function plot_windstorm_event_count(SESSI, varargin)
    % This function takes in the main SESSI table and country-specific tables 
    % as varargin and plots the count of windstorm events (with energy > 0) for each country during boreal winter.
    
    % Filter out windstorm events that have their start time during December, January, or February
    winter_months = [6, 7, 8];
    SESSI.startTime = datetime(SESSI.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    SESSI_winter = SESSI(ismember(month(SESSI.startTime), winter_months),:);
    
    % For each country, count the number of events with energy above 0
    num_countries = length(varargin);
    countries = cell(1, num_countries);
    event_counts = zeros(1, num_countries);
    for i = 1:num_countries
        country_table = varargin{i};
        country_table.startTime = datetime(country_table.startTime, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
        country_winter = country_table(ismember(month(country_table.startTime), winter_months) & country_table.stormEnergy > 0,:);
        event_counts(i) = height(country_winter);
        countries{i} = extractAfter(inputname(i+1), "SE_");  % Extract country code from variable name
    end
    
    % Plot
    figure;
    bar(event_counts);
    ax = gca;
    ax.XTick = 1:num_countries;
    ax.XTickLabel = countries;
    set(gca, 'XTickLabelRotation', 45); % Rotate x-axis labels by 45 degrees for better visibility
    xlabel('Country');
    ylabel('Number of Windstorm Events');
    title('Number of Boreal Winter Windstorm Events (Energy > 0) by Country');
    grid on;
end
