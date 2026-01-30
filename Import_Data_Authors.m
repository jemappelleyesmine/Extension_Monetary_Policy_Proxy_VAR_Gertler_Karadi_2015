% Import_Data_Authors.m
% =========================================================================
% Loads AUTHORS' original datasets:
%   data/karadi_2015_data/VAR_data.csv
%   data/karadi_2015_data/factor_data.csv
%
% Builds DATASET_VAR and DATASET_FACTORS (TSERIES + MAP) exactly as VAR_main expects.
% =========================================================================

fprintf('Loading AUTHORS datasets (Karadi 2015)...\n');

%% Robust paths relative to this script location
thisDir = fileparts(mfilename('fullpath'));

% Expected location: <project_root>/data/karadi_2015_data/
dataDir = fullfile(thisDir, 'data', 'karadi_2015_data');

% If scripts live in /code or similar, allow one-level-up
if ~exist(dataDir, 'dir')
    dataDir = fullfile(thisDir, '..', 'data', 'karadi_2015_data');
end

% If user runs from another working directory, still fail loudly if not found
if ~exist(dataDir, 'dir')
    error('Could not find data folder "data/karadi_2015_data" relative to %s', thisDir);
end

varPath = fullfile(dataDir, 'VAR_data.csv');
facPath = fullfile(dataDir, 'factor_data.csv');

if ~isfile(varPath)
    error('Missing file: %s', varPath);
end
if ~isfile(facPath)
    error('Missing file: %s', facPath);
end

%% Read CSVs
Tvar = readtable(varPath);
Tfac = readtable(facPath);

fprintf('  VAR data: %d observations (%d-%02d to %d-%02d)\n', ...
    height(Tvar), Tvar.year(1), Tvar.month(1), Tvar.year(end), Tvar.month(end));
fprintf('  Factors data: %d observations (%d-%02d to %d-%02d)\n', ...
    height(Tfac), Tfac.year(1), Tfac.month(1), Tfac.year(end), Tfac.month(end));

%% Zero-fill ONLY event-based proxies (if they contain NaNs)
% Do NOT touch ff4_tc: if it is missing outside its availability, keep NaN.
cols_to_zero = {'mp1_tc','ed2_tc','ed3_tc','ed4_tc'};
for c = 1:numel(cols_to_zero)
    v = cols_to_zero{c};
    if ismember(v, Tfac.Properties.VariableNames)
        nan_count = sum(isnan(Tfac.(v)));
        Tfac.(v)(isnan(Tfac.(v))) = 0;
        if nan_count > 0
            fprintf('  Converted %d NaN values to 0 in %s\n', nan_count, upper(v));
        end
    end
end

%% Build dataset structures
DATASET_VAR = struct();
DATASET_FACTORS = struct();

DATASET_VAR.TSERIES     = table2array(Tvar);
DATASET_FACTORS.TSERIES = table2array(Tfac);

%% Build MAPs (case-insensitive by using UPPER)
vnames = upper(Tvar.Properties.VariableNames);
fnames = upper(Tfac.Properties.VariableNames);

DATASET_VAR.MAP     = containers.Map(vnames, 1:numel(vnames));
DATASET_FACTORS.MAP = containers.Map(fnames, 1:numel(fnames));

%% Aliases that GK code expects
% (Authors VAR_data.csv usually already has LCPI/LIP, but add guards anyway.)
if isKey(DATASET_VAR.MAP, 'LOGCPI')
    DATASET_VAR.MAP('LCPI') = DATASET_VAR.MAP('LOGCPI');
end
if isKey(DATASET_VAR.MAP, 'LOGIP')
    DATASET_VAR.MAP('LIP') = DATASET_VAR.MAP('LOGIP');
end

%% Spread aliases (harmless if not used)
if isKey(DATASET_VAR.MAP, 'MORTG_SPREAD_M')
    DATASET_VAR.MAP('MORTG_SPREAD') = DATASET_VAR.MAP('MORTG_SPREAD_M');
end
if isKey(DATASET_VAR.MAP, 'CP3M_SPREAD_M')
    DATASET_VAR.MAP('CP3M_SPREAD') = DATASET_VAR.MAP('CP3M_SPREAD_M');
end

fprintf('AUTHORS data loaded successfully!\n\n');

%% Push to base workspace for VAR_main to use
assignin('base', 'DATASET_VAR', DATASET_VAR);
assignin('base', 'DATASET_FACTORS', DATASET_FACTORS);
