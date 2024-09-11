% S = downscaleContour(S);
EUWS = applyMask(S, EUWS_Raw);


function S = downscaleContour(S)
fields = fieldnames(S);
% Loop over each country
for i = 1:length(fields)
    country = fields{i};

    % Loop over all points in a country
    for j = 1:length(S.(country))
        polygon = [S.(country)(j).X', S.(country)(j).Y'];
        % Remove NaN values
        polygon(any(isnan(polygon), 2), :) = [];
        % Set tolerance
        tolerance = 0.001;
        polygon = reducepoly(polygon, tolerance);
        % Update the polygon data
        S.(country)(j).X = polygon(:, 1)';
        S.(country)(j).Y = polygon(:, 2)';
    end
end

end

function EUWS = applyMask(S, EUWS_Raw)

% Rough filter
EUWS = EUWS_Raw(EUWS_Raw.latitude_rel_vor   <= 70 &...
            EUWS_Raw.latitude_rel_vor       >= 30 &...
            EUWS_Raw.longitude_rel_vor      <= 30 &...
            EUWS_Raw.longitude_rel_vor      >=-30 & ...
            EUWS_Raw.wind_speed_10m         >= 20 , :);


insideCountry = false(height(EUWS), 1);
% Loop over each windstorm
for i = 1:height(EUWS)

    % Loop over each country
    for fieldName = fieldnames(S)'
        country = S.(fieldName{1});

        % Loop over each polygon in the country
        for k = 1:length(country)
            inside = inpolygon(EUWS.longitude_rel_vor(i), EUWS.latitude_rel_vor(i), country(k).X, country(k).Y);

            % If the point is inside this polygon, set insideCountry(i) to true and break the loop
            if inside
                insideCountry(i) = true;
                break;
            end
        end

        % If the point is already found to be inside a country, break the outer loop as well
        if insideCountry(i)
            break;
        end
    end
end
EUWS = EUWS(insideCountry, :);

end





