%unzipAndMove;
S = countryContours;
%EUWS_Raw = importWindstormFiles;

function unzipAndMove
source_dir = 'D:\ERA5_1HR';
% target_dir = 'D:\1day_500km_mslp_wind925_wind10m_nc';
target_dir = 'D:\2day_1000km_mslp_wind925_wind10m_nc';

subdirs = dir(source_dir);
subdirs(~[subdirs.isdir]) = [];  % remove non-directories
subdirs = subdirs(3:end);  % remove '.' and '..'

% Loop through each subdirectory
for i = 1:length(subdirs)
    subdir_path = fullfile(source_dir, subdirs(i).name);
    % gz_files = dir(fullfile(subdir_path, '*.1day-500km_addmslp_addwind925_addwind10m.new.nc.gz'));
    gz_files = dir(fullfile(subdir_path, '*.2day-1000km_addmslp_addwind925_addwind10m.new.nc.gz'));
    source_year = subdirs(i).name(17:31);

    % Loop through each .gz file
    for j = 1:length(gz_files)
        gz_file_path = fullfile(subdir_path, gz_files(j).name);
        gunzip(gz_file_path);
        [~,name,~] = fileparts(gz_file_path);
        unzipped_file_path = fullfile(subdir_path, name);
        % Add source year to the unzipped file name
        new_file_name = [source_year '_' name];
        new_file_path = fullfile(subdir_path, new_file_name);
        movefile(unzipped_file_path, new_file_path);
        % Move the renamed file to the target directory
        movefile(new_file_path, target_dir);
    end
end
end


function S = countryContours
source_dir = 'D:\contour_files\countries';

subdirs = dir(source_dir);
subdirs(~[subdirs.isdir]) = [];  % remove non-directories
subdirs = subdirs(3:end);  % remove '.' and '..'

S = struct();

% Loop through each subdirectory
for i = 1:length(subdirs)
    subdir_path = fullfile(source_dir, subdirs(i).name);
    shapefiles = dir(fullfile(subdir_path, '*_0.shp'));

    % Loop through each shapefile
    for j = 1:length(shapefiles)
        filename = fullfile(subdir_path, shapefiles(j).name);
        data = shaperead(filename);
        country_code = data(1).GID_0;
        S.(['country_' country_code]) = data;
    end
end

% Modify Russia's data to remove parts with 15 < X <= 25
if isfield(S, 'country_RUS')
    for k = 1:length(S.country_RUS)
        idx = S.country_RUS(k).X > 15 & S.country_RUS(k).X <= 25;
        S.country_RUS(k).X = S.country_RUS(k).X(idx);
        S.country_RUS(k).Y = S.country_RUS(k).Y(idx);  % truncate Y based on X
    end
end

end


function EUWS_Raw = importWindstormFiles
path = 'D:\1day_500km_mslp_wind925_wind10m_nc';
filename = 'oct-mar20062007_tr_trs_pos.1day-500km_addmslp_addwind925_addwind10m.new.nc';
fullpath = fullfile(path, filename);

%info = ncinfo(fullpath); 
%ncdisp(fullpath);

infofile = table(ncread(fullpath, 'TRACK_ID'),    ncread(fullpath, 'FIRST_PT'),   ncread(fullpath, 'NUM_PTS'));

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
EUWS_Raw.Windstorm_Name = zeros(height(EUWS_Raw), 1);
windstorm_name = 0;
for i = 2:height(EUWS_Raw)
    if EUWS_Raw.index(i) < EUWS_Raw.index(i-1)
        windstorm_name = windstorm_name + 1;
    end
    EUWS_Raw.Windstorm_Name(i) = windstorm_name;
end

end


