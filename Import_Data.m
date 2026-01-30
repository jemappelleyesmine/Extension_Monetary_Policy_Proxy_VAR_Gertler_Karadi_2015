% Import_Data.m
% Loads extended datasets and builds TSERIES + MAP structures expected by VAR_main.
% Improvements:
%   - Do NOT zero-fill ff4_tc (NaN means "unavailable instrument", not "no meeting")
%   - Only zero-fill event-based proxies (mp1_tc, ed2_tc, ed3_tc, ed4_tc)
%   - Enforce FF4 availability window (1991:01–2012:06) defensively
%   - Add aliases for spread variable names expected by GK scripts
%   - Keep existing aliases (LOGCPI->LCPI, LOGIP->LIP)
%
fprintf('Loading extended datasets...\n');

%% ------------------------------------------------------------------------
% Robust paths (relative to Import_Data.m location, not pwd)
thisDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(thisDir, 'data', 'output_data');

if ~exist(dataDir, 'dir')
    dataDir = fullfile(thisDir, '..', 'data', 'output_data');
end
if ~exist(dataDir, 'dir')
    dataDir = fullfile(thisDir, 'output_data');
end

varPath = fullfile(dataDir, 'VAR_data_extended.csv');
facPath = fullfile(dataDir, 'Factors_data_extended.csv');

if ~isfile(varPath)
    error('Missing file: %s', varPath);
end
if ~isfile(facPath)
    error('Missing file: %s', facPath);
end

%% ------------------------------------------------------------------------
% Read CSV files
Tvar = readtable(varPath);
Tfac = readtable(facPath);

fprintf('  VAR data: %d observations (%d-%02d to %d-%02d)\n', ...
    height(Tvar), Tvar.year(1), Tvar.month(1), Tvar.year(end), Tvar.month(end));
fprintf('  Factors data: %d observations (%d-%02d to %d-%02d)\n', ...
    height(Tfac), Tfac.year(1), Tfac.month(1), Tfac.year(end), Tfac.month(end));

%% ------------------------------------------------------------------------
% Zero-fill ONLY event-based proxies
% Months without FOMC meetings => 0 surprises for ED2/ED3/ED4 (and mp1 if present).
% DO NOT zero-fill ff4_tc: NaN means instrument not available outside its window.
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

%% ------------------------------------------------------------------------
% Defensive enforcement of FF4 window (even if upstream changes)
if ismember('ff4_tc', Tfac.Properties.VariableNames)
    d = datetime(Tfac.year, Tfac.month, 1);
    ff4_start = datetime(1991,1,1);
    ff4_end   = datetime(2012,6,1);

    % Keep FF4 only inside the availability window; outside => NaN
    out_count = sum(~isnan(Tfac.ff4_tc) & (d < ff4_start | d > ff4_end));
    Tfac.ff4_tc(d < ff4_start | d > ff4_end) = NaN;

    if out_count > 0
        fprintf('  Set %d FF4 values to NaN outside 1991:01–2012:06 (defensive window)\n', out_count);
    end
else
    fprintf('  Note: ff4_tc not found in Factors dataset.\n');
end

%% ------------------------------------------------------------------------
% Build TSERIES matrices exactly as stored in CSVs
DATASET_VAR = struct();
DATASET_FACTORS = struct();

DATASET_VAR.TSERIES     = table2array(Tvar);
DATASET_FACTORS.TSERIES = table2array(Tfac);

%% ------------------------------------------------------------------------
% Build MAPs (case-insensitive by using UPPER)
vnames = upper(Tvar.Properties.VariableNames);
fnames = upper(Tfac.Properties.VariableNames);

DATASET_VAR.MAP     = containers.Map(vnames, 1:numel(vnames));
DATASET_FACTORS.MAP = containers.Map(fnames, 1:numel(fnames));

%% ------------------------------------------------------------------------
% Aliases that GK code expects (logcpi → LCPI, logip → LIP)
if isKey(DATASET_VAR.MAP, 'LOGCPI')
    DATASET_VAR.MAP('LCPI') = DATASET_VAR.MAP('LOGCPI');
end
if isKey(DATASET_VAR.MAP, 'LOGIP')
    DATASET_VAR.MAP('LIP') = DATASET_VAR.MAP('LOGIP');
end

%% ------------------------------------------------------------------------
% Aliases for spread names expected by GK scripts
% Your build script outputs: mortg_spread_m, cp3m_spread_m
% GK scripts often use:      MORTG_SPREAD, CP3M_SPREAD
if isKey(DATASET_VAR.MAP, 'MORTG_SPREAD_M')
    DATASET_VAR.MAP('MORTG_SPREAD') = DATASET_VAR.MAP('MORTG_SPREAD_M');
end
if isKey(DATASET_VAR.MAP, 'CP3M_SPREAD_M')
    DATASET_VAR.MAP('CP3M_SPREAD') = DATASET_VAR.MAP('CP3M_SPREAD_M');
end

%% ------------------------------------------------------------------------
% Optional: quick diagnostics that help catch silent issues early
if isKey(DATASET_FACTORS.MAP, 'FF4_TC')
    ff4_col = DATASET_FACTORS.MAP('FF4_TC');
    ff4_vec = DATASET_FACTORS.TSERIES(:, ff4_col);
    fprintf('  FF4 availability: non-missing months = %d (of %d)\n', ...
        sum(~isnan(ff4_vec)), numel(ff4_vec));
end

fprintf('Data loaded successfully!\n');
fprintf('  Replication period (intended): 1979:07 - 2012:06\n');
fprintf('  Extended period available:     1979:07 - %d:%02d\n\n', Tvar.year(end), Tvar.month(end));

% Push to base workspace for VAR_main to use
assignin('base', 'DATASET_VAR', DATASET_VAR);
assignin('base', 'DATASET_FACTORS', DATASET_FACTORS);
