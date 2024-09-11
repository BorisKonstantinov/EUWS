% SESSI = matchEU(SE, SSIFULL);
% plotStormEnergyVsEU(SESSI);
% SESSI = addNAOIndex(SESSI, naoDpM);

% SE_AUT = addNAOIndex(SE_AUT, naoDpM);
% SE_CHE = addNAOIndex(SE_CHE, naoDpM);
% SE_CZE = addNAOIndex(SE_CZE, naoDpM);
% SE_DEU = addNAOIndex(SE_DEU, naoDpM);
% SE_DNK = addNAOIndex(SE_DNK, naoDpM);
% SE_EST = addNAOIndex(SE_EST, naoDpM);
% SE_FIN = addNAOIndex(SE_FIN, naoDpM);
% SE_FRA = addNAOIndex(SE_FRA, naoDpM);
% SE_GBR = addNAOIndex(SE_GBR, naoDpM);
% SE_HUN = addNAOIndex(SE_HUN, naoDpM);
% SE_IRL = addNAOIndex(SE_IRL, naoDpM);
% SE_LTU = addNAOIndex(SE_LTU, naoDpM);
% SE_LVA = addNAOIndex(SE_LVA, naoDpM);
% SE_NLD = addNAOIndex(SE_NLD, naoDpM);
% SE_NOR = addNAOIndex(SE_NOR, naoDpM);
% SE_POL = addNAOIndex(SE_POL, naoDpM);
% SE_SVK = addNAOIndex(SE_SVK, naoDpM);
% SE_SWE = addNAOIndex(SE_SWE, naoDpM);



function SE = matchEU(SE, SSIFULL)

% Convert EventID to datetime format
SSIFULL.EventID = datetime(string(SSIFULL.EventID), 'InputFormat', 'yyyyMMdd');

% Initialize the new column in SE
SE.EU = NaN(height(SE), 1);

% Iterate over rows in SE
for i = 1:height(SE)
    % Get start and end time for current row
    startTime = SE.startTime(i);
    endTime = SE.endTime(i);

    % Find rows in SSIFULL where EventID is within startTime and endTime
    rows = SSIFULL.EventID >= startTime & SSIFULL.EventID <= endTime;

    % If any matching rows are found
    if any(rows)
        % Take the mean of corresponding EU values and add to SE
        SE.EU(i) = sum(SSIFULL.EU(rows));
    end
end
end

function plotStormEnergyVsEU(SESSI)
    % Ensure there are no rows with NaN values
    SESSI = rmmissing(SESSI);

    % Normalize the stormEnergy and EU columns using min-max normalization
    SESSI.stormEnergy = (SESSI.stormEnergy - min(SESSI.stormEnergy)) / (max(SESSI.stormEnergy) - min(SESSI.stormEnergy));
    SESSI.EU = (SESSI.EU - min(SESSI.EU)) / (max(SESSI.EU) - min(SESSI.EU));

    % Create a new figure
    figure;

    % Create a scatter plot of normalized stormEnergy vs EU
    scatter(SESSI.stormEnergy, SESSI.EU, 'filled');

    % Fit a linear regression model and get the predicted values
    lm = fitlm(SESSI, 'EU ~ stormEnergy');
    hold on;
    plot(SESSI.stormEnergy, lm.Fitted, 'r');
    hold off;

    % Set the title and labels
    title('Normalized Storm Energy vs EU');
    xlabel('Normalized Storm Energy');
    ylabel('Normalized EU');

    % Add a grid for easier visualization
    grid on;

    % Get the coefficients from the fit
    coeff = lm.Coefficients.Estimate;
    
    % Create the equation string
    equation = ['y = ', num2str(coeff(2)), 'x + ', num2str(coeff(1))];
    
    % Get the R-squared value
    R2 = lm.Rsquared.Ordinary;
    
    % Create the R-squared string
    Rsqrd = ['R^2 = ', num2str(R2)];

    % Add the equation and the R-squared value to the plot
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', {equation, Rsqrd}, 'FitBoxToText', 'on');
end

function SESSI = addNAOIndex(SESSI, naoDpM)
    % Initialize new column in SE
    SESSI.NAOIndex = NaN(height(SESSI), 1);

    % For each row in SE, find matching index_m value(s) in naoDpM
    for i = 1:height(SESSI)
        % Get start and end times for this row
        startTime = SESSI.startTime(i);
        endTime = SESSI.endTime(i);

        % Find rows in naoDpM where the date is within start and end times
        rows = (naoDpM.date >= startTime) & (naoDpM.date <= endTime);

        % If any matching rows were found, take the average index_m value
        if any(rows)
            SESSI.NAOIndex(i) = mean(naoDpM.index_m(rows));
        end
    end
end


