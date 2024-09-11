clc;
path = 'D:\1day_500km_mslp_wind925_wind10m_nc';
filename = 'apr-sep2020_DET_tr_trs_pos.1day-500km_addmslp_addwind925_addwind10m.new.nc';
fullpath = fullfile(path, filename);
%info = ncinfo(fullpath); 
%ncdisp(fullpath);

infofile = table(ncread(fullpath, 'TRACK_ID'),    ncread(fullpath, 'FIRST_PT'),   ncread(fullpath, 'NUM_PTS'));

% Create EUWS table
EUWS =  table(ncread(fullpath, 'index'),       ncread(fullpath, 'time'),...
              ncread(fullpath, 'longitude'),   ncread(fullpath, 'latitude'),   ncread(fullpath, 'relative_vorticity'),...
              ncread(fullpath, 'longitude_1'), ncread(fullpath, 'latitude_1'), ncread(fullpath, 'air_pressure_at_sea_level'),...
              ncread(fullpath, 'longitude_2'), ncread(fullpath, 'latitude_2'), ncread(fullpath, 'wind_speed_925'),...
              ncread(fullpath, 'longitude_3'), ncread(fullpath, 'latitude_3'), ncread(fullpath, 'wind_speed_10m'));

EUWS.Properties.VariableNames = {'index',               'time',...
                                 'longitude_rel_vor',   'latitude_rel_vor', 'relative_vorticity',...
                                 'longitude_apsl',      'latitude_apsl',    'apsl',...
                                 'longitude_ws925',     'latitude_ws925',   'wind_speed_925',...
                                 'longitude_ws10m',     'latitude_ws10m',   'wind_speed_10m'};
origin = datetime(1979, 1, 1, 0, 0, 0);  % Origin time
EUWS.time = origin + hours(EUWS.time);  % Convert hours since origin to datetime
EUWS.date = dateshift(EUWS.time, 'start', 'day');



% Filter
EUWS = EUWS(EUWS.latitude_rel_vor   <= 70 &...
            EUWS.latitude_rel_vor   >= 30 &...
            EUWS.longitude_rel_vor  <= 30 &...
            EUWS.longitude_rel_vor  >=-30 & ...
            EUWS.wind_speed_10m     >= 20 , :);

% Initialize a logical array to store results of inpolygon
insideCountry = false(height(EUWS), 1);

% Loop over each windstorm
for i = 1:height(EUWS)
    % Loop over each country
    for fieldName = fieldnames(S)'
        country = S.(fieldName{1});
        % Loop over each polygon in the country
        for k = 1:length(country)
            % Test if the point is inside the polygon
            inside = inpolygon(EUWS.longitude_ws10m(i), EUWS.latitude_ws10m(i), country(k).X, country(k).Y);
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

% Filter the EUWS table to include only points that are inside a country
EUWS = EUWS(insideCountry, :);



