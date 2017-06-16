%% LEMS_Complete_Averaging.m
% Nipun Gunawardena
% fileConcatenate.py concatenates all LEMS csv files into LEMSx_Latest.CSV
% using python. This script then retrieves the 'Latest' files, averages
% them, and saves it as a .mat file

clear all, close all, clc


%% Get files
disp('Seraching for files...')
latestFiles = {};
fileList = dir();
for i = 1:length(fileList)
    currentName = fileList(i).name;
    if strfind(currentName,'Latest') & (strcmp(currentName(end-2:end),'CSV')) & (~strcmp(currentName(1),'.'))
        latestFiles{end+1} = currentName;
        fprintf('Found %s\n', currentName);
    end
end
numFiles = length(latestFiles);
lemsNames = cell(size(latestFiles));
for i = 1:numFiles
    lemsNames{i} = latestFiles{i}(1:5);
end



%% Get averaging period from user
avgPeriod = input('\nPlease enter averaging period in minutes (recommended value is 5): ');



%% Load Data
disp('Reading files into memory, may take some time...')
lemsRawData = cell(numFiles, 1);
for i = 1:numFiles
    fprintf('Reading %s\n', lemsNames{i});
    lemsRawData{i} = importfile(latestFiles{i});
end
disp('Finished reading files...')



%% Fix LEMS J soil
% LEMS J has switched soil sensors until 2017/01/13/15/54/10 (Index 259295)
% Switch them back
switchIdx = 259295;

tempSoil = lemsRawData{10}.Upper_Soil_Mois(1:switchIdx);
lemsRawData{10}.Upper_Soil_Mois(1:switchIdx) = lemsRawData{10}.Lower_Soil_Mois(1:switchIdx);
lemsRawData{10}.Lower_Soil_Mois(1:switchIdx) = tempSoil;

tempSoil = lemsRawData{10}.Upper_Soil_Temp(1:switchIdx);
lemsRawData{10}.Upper_Soil_Temp(1:switchIdx) = lemsRawData{10}.Lower_Soil_Temp(1:switchIdx);
lemsRawData{10}.Lower_Soil_Temp(1:switchIdx) = tempSoil;



%% Save Raw Data
disp('Saving raw data...')
save('LEMS_Raw_Data.mat', 'lemsRawData');



%% Calculate Wind Components
disp('Calculating wind components...')
for i = 1:numFiles
    cartDeg = 270 - lemsRawData{i}.Wind_Dir;
    cartRad = cartDeg*(pi/180);
    [windU, windV] = pol2cart(cartRad, lemsRawData{i}.Wind_Spd);
    lemsRawData{i}.windU = windU;
    lemsRawData{i}.windV = windV;
end



%% Calculate Virtual Potential Temperature
disp('Calculating virtual potential temperature @ 2m...')
eso = 6.11;         % hPa - Reference saturation vapor pressure
lv = 2.5e6;         % (J/kg)   - Latent heat of vaporization of water
Rv = 461.5;         % (J*K)/kg - Gas constant for water vapor
T0 = 273.15;        % K - Reference temperature
gamma = 0.286;      % R/Cp, where Cp is specific heat at sea level & 300K
for i = 1:numFiles
    T = lemsRawData{i}.SHT_Amb_C + 273.15;  % Convert to K
    P = lemsRawData{i}.Pressure / 100;      % Convert to hPa
    RH = lemsRawData{i}.SHT_Hum_Pct;        % Percentage, no conversion
    es = eso*exp( (lv/Rv) * (1/T0 - 1./T) );
    e = (RH / 100) .* es;
    q = 0.622 * (e./P);
    thetaV = T .* (1.0 + 0.61*q) .* (1000./P).^gamma;
    lemsRawData{i}.thetaV = thetaV;
end



%% Calculate Virtual Potential Temperature for Surface
disp('Calculating virtual potential temperature @ surface...')
pHeights = [0.83, 0.95, 0.98, 0.95, 1.03, 0.97, 0.82, 0.93, 0.91, 0.85, 0.85, 0.97];  % Height of each pressure sensor in m
groundAlt = [332, 347, 397, 328, 293, 383, 325, 276, 385, 262, 317, 368]; % Altitude at each LEMS in m
pSea = 101325;          % Pa
Rs = 8.3144598;         % J/(mol*K)
g0 = 9.80665;           % m/s^2
M = 0.0289644;          % kg/mol
baroForm = @(h, Pb, T) Pb.*exp(-(g0*M*h)./(Rs.*T));

eso = 6.11;         % hPa - Reference saturation vapor pressure
lv = 2.5e6;         % (J/kg)   - Latent heat of vaporization of water
Rv = 461.5;         % (J*K)/kg - Gas constant for water vapor
T0 = 273.15;        % K - Reference temperature
gamma = 0.286;      % R/Cp, where Cp is specific heat at sea level & 300K
for i = 1:numFiles
    T = lemsRawData{i}.MLX_IR_C + 273.15;   % Convert to K
    Pnc = lemsRawData{i}.Pressure;          % Uncorrected pressure
    h = pHeights(i);                        % Height of pressure sensor
    P = baroForm(-h, Pnc, T);               % Correct pressure to ground height
    P = P / 100;                            % Convert to hPa
    RH = lemsRawData{i}.SHT_Hum_Pct;        % Percentage, no conversion
    es = eso*exp( (lv/Rv) * (1/T0 - 1./T) );
    e = (RH / 100) .* es;
    q = 0.622 * (e./P);
    thetaVSurface = T .* (1.0 + 0.61*q) .* (1000./P).^gamma;
    lemsRawData{i}.thetaVSurface = thetaVSurface;
end



%% Find minimum and maximum dates
minDate = 1E9;
maxDate = 0;
for i = 1:numFiles
    lemsRawData{i}.datenum = datenum(lemsRawData{i}.Year, lemsRawData{i}.Month, lemsRawData{i}.Date, lemsRawData{i}.Hour, lemsRawData{i}.Minute, lemsRawData{i}.Second);
    miD = min(lemsRawData{i}.datenum);
    maD = max(lemsRawData{i}.datenum);
    if miD < minDate
        minDate = miD;
    end
    if maD > maxDate
        maxDate = maD;
    end
end



%% Create unified average date vector
minDate = floor(minDate);
maxDate = ceil(maxDate);
dates = zeros((maxDate - minDate)/datenum(0,0,0,0,avgPeriod,0), 1);
dates(1) = minDate;
for i = 2:length(dates)
    dates(i) = addtodate(dates(i-1), avgPeriod, 'minute');
    if (floor(dates(i)) - floor(dates(i-1))) ~= 0
        dates(i) = floor(dates(i));  % Helps prevent rounding errors
    end
end
avgLen = length(dates);

lemsMeasGood = nan(avgLen, numFiles);   % Simple array to indicate whether the LEMS has an average at a given time stamp



%% Create UTC dates vector
utcDates = dates - datenum([0,0,0,1,0,0]);



%% Create cell array of tables for 5 min average data
[~, rawNumCols] = size(lemsRawData{1});
lemsAvgData = cell(numFiles, 1);    % Each cell will hold a table
for i = 1:numFiles
    lemsAvgData{i} = table();
    lemsAvgData{i}.dates = dates;
    lemsAvgData{i}.utcDates = utcDates;
    for col = 7:(rawNumCols-1)  % The last column of lemsRawData is the datenum col. If variables that need to be averaged are to be added, do it before the 'Find minimum and maximum dates' section, or bugs will happen
        lemsAvgData{i}.(lemsRawData{i}.Properties.VariableNames{col}) = nan(avgLen, 1);
    end
end
[~, avgNumCols] = size(lemsAvgData{1});



%% Actually perform averages
disp(' ');  % Formatting
disp('Beginnging Averaging...');
tTotal = 0;
for i = 1:numFiles
    fprintf('Beginning Averaging %s...\n', lemsNames{i});
    lemsAvgData{i}.Properties.Description = lemsNames{i};
    tStart = tic;
    cLems = lemsRawData{i};
    lowBound = find(dates <= cLems.datenum(1), 1, 'last')+1;
    highBound = find(dates >= cLems.datenum(end), 1, 'first');
    startTime = dates(lowBound - 1);
    endTime = dates(lowBound);
    for j = lowBound:highBound
        valid = (cLems.datenum > startTime) & (cLems.datenum <= endTime);
        for k = 3:avgNumCols
            lemsAvgData{i}(j,k) = {nanmean(table2array(cLems(valid, k+4)))};   % +4 is to get to right column. cLems has extra date variables lemsAvgData doesn't
        end
        if ~isnan(lemsAvgData{i}{j,k})
            lemsMeasGood(j,i) = i;
        end
        startTime = endTime;
        endTime = dates(j+1);
    end
    lemsAvgData{i}.Properties.VariableUnits = {'Matlab Datenum','Matlab Datenum','Volts','Celsius','Celsius','Celsius','meters^3 H20/meters^3 soil','Celsius','meters^3 H20/meters^3 soil','Pascals','Celsius','Degrees from north','meters/second','Watts/meter^2','Celsius','% Relative Humidity','meters/second','meters/second', 'Kelvin', 'Kelvin'};
    lemsAvgData{i}.Properties.VariableDescriptions = {'Measurement timestamp Local','Measurement timestamp UTC','LEMS battery/charger level','MLX90614 surface temperature','MLX90614 ambient temperature','5TM shallow soil temperature','5TM shallow soil VWC','5TM deep soil temperature','5tm deep soil VWC','BMP280 barometric pressure','BMP280 box temperature','Davis wind direction','Davis wind speed','LI200-R incoming shortwave radiation','SHT31 ambient temperature','SHT31 ambient relative humidity','Davis wind U component','Davis wind V component', 'Virtual Potential Temperature @ 2m', 'Virtual Potential Surface Temperature'};
    fprintf('Finished Averaging %s.\n', lemsNames{i});
    tElapsed = toc(tStart);
    fprintf('Elapsed averaging time for this LEMS is %f seconds\n', tElapsed);
    tTotal = tTotal + tElapsed;
    disp(' ');
end
fprintf('The total averaging time is %f hours\n', tTotal/3600);



%% Remove bad values
for i = 1:numFiles
    % Upper Soil Moisture
    valid = lemsAvgData{i}.Upper_Soil_Mois < 0;
    lemsAvgData{i}.Upper_Soil_Mois(valid) = NaN;
    
    % Lower Soil Moisture
    valid = lemsAvgData{i}.Lower_Soil_Mois < 0;
    lemsAvgData{i}.Lower_Soil_Mois(valid) = NaN;
    
    % Upper Soil Temperature
    valid = lemsAvgData{i}.Upper_Soil_Temp < -200.0;
    lemsAvgData{i}.Upper_Soil_Temp(valid) = NaN;
    
    % Lower Soil Temperature
    valid = lemsAvgData{i}.Lower_Soil_Temp < -200.0;
    lemsAvgData{i}.Lower_Soil_Temp(valid) = NaN;
end



%% Save In Case Want to Load
save('LEMS_Avg_Latest.mat', 'lemsAvgData', 'lemsMeasGood', 'dates', 'minDate', 'maxDate', 'lemsNames', 'numFiles');



%% Signify averaging is done
% load('train.mat');
% soundsc(y);