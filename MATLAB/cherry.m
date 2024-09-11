% mesh = gridMaker(lon,lat);
% maskedMatrix = createMask(S, mesh);

% maskedMatrixSWE = createMask(S_SWE, mesh);
% logicMapOfSWE = flipud(maskedMatrixSWE);

% Run below when wanting to reset SE. This means EnergyCalculator can be
% run in segments of 200 for each of the 7 parts.
% SE = EUWS_DateRange;
% SE.stormEnergy = zeros(height(SE), 1);

% SE = EnergyCalculator(SE, logicMapOfEU);

% For a country
% SE = EUWS_DateRange;
% SE.stormEnergy = zeros(height(SE), 1);

% SE_AUT = EnergyCalculator(SE_FIC, logicMapOfAUT);
% SE_CHE = EnergyCalculator(SE_FIC, logicMapOfCHE);
% SE_CZE = EnergyCalculator(SE_FIC, logicMapOfCZE);
% SE_DEU = EnergyCalculator(SE_FIC, logicMapOfDEU);
% SE_DNK = EnergyCalculator(SE_FIC, logicMapOfDNK);
% SE_EST = EnergyCalculator(SE_FIC, logicMapOfEST);
% SE_FIN = EnergyCalculator(SE_FIC, logicMapOfFIN);
% SE_FRA = EnergyCalculator(SE_FIC, logicMapOfFRA);
% SE_GBR = EnergyCalculator(SE_FIC, logicMapOfGBR);
% SE_HUN = EnergyCalculator(SE_FIC, logicMapOfHUN);
% SE_IRL = EnergyCalculator(SE_FIC, logicMapOfIRL);
% SE_LTU = EnergyCalculator(SE_FIC, logicMapOfLTU);
% SE_LVA = EnergyCalculator(SE_FIC, logicMapOfLVA);
% SE_NLD = EnergyCalculator(SE_FIC, logicMapOfNLD);
% SE_NOR = EnergyCalculator(SE_FIC, logicMapOfNOR);
% SE_POL = EnergyCalculator(SE_FIC, logicMapOfPOL);
% SE_SVK = EnergyCalculator(SE_FIC, logicMapOfSVK);
% SE_SWE = EnergyCalculator(SE_FIC, logicMapOfSWE);

energyHeatmapNAOPos = EnergyCalculatorHeatmap(SESSI)



function SE = EnergyCalculator(SE, logicMapOfEU)


long = ncread('D:\ERA5-Land\1U_ERA5LAND.nc', 'lon');
latd = ncread('D:\ERA5-Land\1U_ERA5LAND.nc', 'lat');
latd = flip(latd);



% FOR FILESET 1
time = ncread('D:\ERA5-Land\1U_ERA5LAND.nc', 'time');
origin = datetime(1950, 3, 1, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 1:200
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\1U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\1V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end



% FOR FILESET 2
time = ncread('D:\ERA5-Land\2U_ERA5LAND.nc', 'time');
origin = datetime(1961, 11, 2, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 201:400
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\2U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\2V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end



% FOR FILESET 3
time = ncread('D:\ERA5-Land\3U_ERA5LAND.nc', 'time');
origin = datetime(1972, 3, 25, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 401:600
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\3U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\3V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end



% FOR FILESET 4
time = ncread('D:\ERA5-Land\4U_ERA5LAND.nc', 'time');
origin = datetime(1982, 4, 1, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 601:800
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\4U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\4V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end



% FOR FILESET 5
time = ncread('D:\ERA5-Land\5U_ERA5LAND.nc', 'time');
origin = datetime(1992, 10, 23, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 801:1000
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\5U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\5V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end



% FOR FILESET 6
time = ncread('D:\ERA5-Land\6U_ERA5LAND.nc', 'time');
origin = datetime(2004, 5, 2, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 1001:1200
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\6U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\6V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end



% FOR FILESET 7
time = ncread('D:\ERA5-Land\7U_ERA5LAND.nc', 'time');
origin = datetime(2015, 9, 1, 0, 0, 0);   % Origin time for 1
time = origin + hours(time);  % Convert hours since origin to datetime

for DRI = 1201:1322
    DRI
    startIndex = find(time == SE.startTime(DRI));
    endIndex = startIndex + hours(SE.duration(DRI));

    for tP = startIndex:endIndex
        chunkU = ncread('D:\ERA5-Land\7U_ERA5LAND.nc', 'u10', [1, 1, tP], [471, 301, 1]);
        chunkU = rot90(chunkU, 1);
        chunkV = ncread('D:\ERA5-Land\7V_ERA5LAND.nc', 'v10', [1, 1, tP], [471, 301, 1]);
        chunkV = rot90(chunkV, 1);

        for j = 1:301
            for i = 1:471
                if logicMapOfEU(j, i) == 1
                    ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                    if ws >= 15
                        Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                        Energy = 1.2 * Length*10 * ws^3 * 60*60;
                        SE.stormEnergy(DRI) = SE.stormEnergy(DRI) + Energy;
                    end
                end
            end
        end
    end
end


end


function mesh = gridMaker(lon,lat)
% assuming lat and lon are your arrays
[lonGrid, latGrid] = meshgrid(lon, lat);

% dimensions
[nLat, nLon] = size(latGrid);

% initialize cell array
geoCoordsCellArray = cell(nLat, nLon);

% fill cell array
for i = 1:nLat
    for j = 1:nLon
        geoCoordsCellArray{i, j} = [latGrid(i, j), lonGrid(i, j)];
    end
end
mesh = geoCoordsCellArray;

end


function maskedMatrix = createMask(S, coordinateMatrix)
insideCountry = false(301, 471);

for j = 1:height(coordinateMatrix)
    for i = 1:width(coordinateMatrix)
        cell = coordinateMatrix{j, i};
        cellX = cell(2);
        cellY = cell(1);
        for fieldName = fieldnames(S)'
            country = S.(fieldName{1});
            for k = 1:length(country)
                inside = inpolygon(cellX, cellY, country(k).X, country(k).Y);

                if inside
                    insideCountry(j, i) = true;
                    break;
                end
            end

            if insideCountry(j, i)
                break;
            end
        end
    end
end

maskedMatrix = insideCountry;

end


function Length = Haversine(lon, lat, u, v)
% Add 90 to theta to get a line normal to the wind
theta = atan2d(v, u) + 90;

if ( (theta <= 45 && theta >= -45) || theta >= 135 || theta <= -135)
    theta = 90 - abs(mod(theta+90, 360) -180);
    lon1 = lon + 0.05*(sign(u)+ (u==0));
    lat1 = lat + 0.05*tand(theta);
    lon2 = lon - 0.05*(sign(u)+ (u==0));
    lat2 = lat - 0.05*tand(theta);
else
    theta = 45 - abs(mod(abs(theta)+90, 270) -135);
    lon1 = lon + 0.05*tand  (theta);
    lat1 = lat + 0.05*(sign(v)+ (v==0));
    lon2 = lon - 0.05*tand(theta);
    lat2 = lat - 0.05*(sign(v)+ (v==0));
end
% a = sin²(Δφ/2) + cos φ1 ⋅ cos φ2 ⋅ sin²(Δλ/2)
% where	φ is latitude, λ is longitude in radians
a = sind(abs(lat1 - lat2)/2)^2 + cosd(lat1)*cosd(lat2)*sind(abs(lon1 - lon2)/2)^2;
% c = 2 ⋅ atan2( √a, √(1−a) )
c = 2*atan2(sqrt(a), sqrt(1-a));
% d = R ⋅ c
% where R is the radius of the Earth = 6371km
Length = 6371000 * c;
end


function energyHeatmap = EnergyCalculatorHeatmap(SESSI)

long = ncread('D:\ERA5-Land\1U_ERA5LAND.nc', 'lon');
latd = ncread('D:\ERA5-Land\1U_ERA5LAND.nc', 'lat');
latd = flip(latd);

energyHeatmap = zeros(301,471);

origins = [datetime(1950, 3, 1, 0, 0, 0), datetime(1961, 11, 2, 0, 0, 0), datetime(1972, 3, 25, 0, 0, 0), datetime(1982, 4, 1, 0, 0, 0), datetime(1992, 10, 23, 0, 0, 0), datetime(2004, 5, 2, 0, 0, 0), datetime(2015, 9, 1, 0, 0, 0)];
fileSets = 1:7;
ranges = {[1,200], [201,400], [401,600], [601,800], [801,1000], [1001,1200], [1201,1318]};

for fs = fileSets
    time = ncread(['D:\ERA5-Land\', num2str(fs), 'U_ERA5LAND.nc'], 'time');
    origin = origins(fs);
    time = origin + hours(time);
    
    for DRI = ranges{fs}(1):ranges{fs}(2)
        DRI
        if SESSI.NAOIndex(DRI) < -0.5
        %if (SESSI.NAOIndex(DRI) >= -0.5 && SESSI.NAOIndex(DRI) <= 0.5)
            startIndex = find(time == SESSI.startTime(DRI));
            endIndex = startIndex + hours(SESSI.duration(DRI));
            
            for tP = startIndex:endIndex
                chunkU = ncread(['D:\ERA5-Land\', num2str(fs), 'U_ERA5LAND.nc'], 'u10', [1, 1, tP], [471, 301, 1]);
                chunkU = rot90(chunkU, 1);
                chunkV = ncread(['D:\ERA5-Land\', num2str(fs), 'V_ERA5LAND.nc'], 'v10', [1, 1, tP], [471, 301, 1]);
                chunkV = rot90(chunkV, 1);
                
                for j = 1:301
                    for i = 1:471
                        ws = sqrt(chunkU(j, i)^2 + chunkV(j, i)^2);
                        if ws >= 15
                            Length = Haversine(long(i), latd(j), chunkU(j, i), chunkV(j, i));
                            Energy = 1.2 * Length*10 * ws^3 * 60*60;
                            energyHeatmap(j, i) = energyHeatmap(j, i) + Energy;
                        end
                    end
                end
            end
        end
    end
end

end
