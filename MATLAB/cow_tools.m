% Step 1: Import raw data from Nov to Sep to EUWS_Raw
% Step 2: Filter data using filters in geomask
% Step 3: Append filtered results to EUWS
% Step 4: Repeat Step 1->3 for raw data from Oct to Mar
% Step 5: Repeat Step 4 for next year


EUWS = wsfromto(1950, 2020, S);


function EUWS = wsfromto(from, to, S)
EUWS = []; % Init. empty table
ws_id = 0; % Init. windstorm name
for year = from:to
    ws_id
    for isOctToMar = [false true]
        % Step 1: Import raw data
        filename = getFilename(year, isOctToMar)
        [EUWS_Raw, ws_id] = importWindstormFiles(filename, ws_id);

        % Step 2: Filter data
        EUWS_Filtered = applyMask(S, EUWS_Raw);

        % Step 3: Append filtered results to EUWS
        EUWS = [EUWS; EUWS_Filtered];
    end
end

end


function filename = getFilename(year, isOctToMar)
if isOctToMar
    nextYear = year + 1;
    filename = ['oct-mar', num2str(year), num2str(nextYear), '_tr_trs_pos.1day-500km_addmslp_addwind925_addwind10m.new.nc'];
else
    filename = ['apr-sep', num2str(year), '_DET_tr_trs_pos.1day-500km_addmslp_addwind925_addwind10m.new.nc'];
end
end


function [EUWS_Raw, ws_id] = importWindstormFiles(filename, ws_id)
path = 'D:\1day_500km_mslp_wind925_wind10m_nc';
fullpath = fullfile(path, filename);

% Create EUWS table
EUWS_Raw =  table(ncread(fullpath, 'index'),       ncread(fullpath, 'time'),...
              ncread(fullpath, 'longitude'),   ncread(fullpath, 'latitude'),   ncread(fullpath, 'relative_vorticity'),...
              ncread(fullpath, 'longitude_1'), ncread(fullpath, 'latitude_1'), ncread(fullpath, 'air_pressure_at_sea_level'),...
              ncread(fullpath, 'longitude_2'), ncread(fullpath, 'latitude_2'), ncread(fullpath, 'wind_speed_925'),...
              ncread(fullpath, 'longitude_3'), ncread(fullpath, 'latitude_3'), ncread(fullpath, 'wind_speed_10m'));

EUWS_Raw.Properties.VariableNames = {'index',               'time',...
                                 'longitude_rel_vor',   'latitude_rel_vor', 'relative_vorticity',...
                                 'longitude_apsl',      'latitude_apsl',    'apsl',...
                                 'longitude_ws925',     'latitude_ws925',   'wind_speed_925',...
                                 'longitude_ws10m',     'latitude_ws10m',   'wind_speed_10m'};
origin = datetime(1979, 1, 1, 0, 0, 0);  % Origin time
EUWS_Raw.time = origin + hours(EUWS_Raw.time);  % Convert hours since origin to datetime
EUWS_Raw.date = dateshift(EUWS_Raw.time, 'start', 'day');

% Modify longitude fields
EUWS_Raw.longitude_rel_vor(EUWS_Raw.longitude_rel_vor > 180) = EUWS_Raw.longitude_rel_vor(EUWS_Raw.longitude_rel_vor > 180) - 360;
EUWS_Raw.longitude_apsl(EUWS_Raw.longitude_apsl > 180) = EUWS_Raw.longitude_apsl(EUWS_Raw.longitude_apsl > 180) - 360;
EUWS_Raw.longitude_ws925(EUWS_Raw.longitude_ws925 > 180) = EUWS_Raw.longitude_ws925(EUWS_Raw.longitude_ws925 > 180) - 360;
EUWS_Raw.longitude_ws10m(EUWS_Raw.longitude_ws10m > 180) = EUWS_Raw.longitude_ws10m(EUWS_Raw.longitude_ws10m > 180) - 360;

% Add Windstorm Name
EUWS_Raw.windstorm_name = zeros(height(EUWS_Raw), 1);
for i = 2:height(EUWS_Raw)
    if EUWS_Raw.index(i) < EUWS_Raw.index(i-1)
        ws_id = ws_id + 1;
    end
    EUWS_Raw.windstorm_name(i) = ws_id;
end

end


function EUWS = applyMask(S, EUWS_Raw)

% Rough filter
EUWS = EUWS_Raw(EUWS_Raw.latitude_rel_vor   <= 72 &...
            EUWS_Raw.latitude_rel_vor       >= 42 &...
            EUWS_Raw.longitude_rel_vor      <= 32 &...
            EUWS_Raw.longitude_rel_vor      >=-12 & ...
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