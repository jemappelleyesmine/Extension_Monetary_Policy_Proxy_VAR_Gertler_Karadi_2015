%% build_extended_datasets.m 
% =========================================================================

clear; clc; close all;

fprintf('========================================================================\n');
fprintf('BUILDING EXTENDED DATASETS\n');
fprintf('========================================================================\n\n');

%% ========================================================================
%  0) PATHS
% =========================================================================
folders.base   = pwd;
folders.data   = fullfile(folders.base, 'data');
folders.raw    = fullfile(folders.data, 'raw_data');
folders.auth15 = fullfile(folders.data, 'karadi_2015_data');
folders.output = fullfile(folders.data, 'output_data');
if ~exist(folders.output,'dir'); mkdir(folders.output); end

% Raw inputs
files.cpi      = fullfile(folders.raw, 'CPIAUCSL.csv');
files.ip       = fullfile(folders.raw, 'INDPRO.csv');
files.gs1      = fullfile(folders.raw, 'GS1.csv');
files.gsw      = fullfile(folders.raw, 'feds200628.csv'); % GSW curve (daily or monthly)
files.gz_m     = fullfile(folders.raw, 'GZ_monthly.csv'); % GZ monthly with ebp (1973-2010ish)
files.ebp_ext  = fullfile(folders.raw, 'ebp_csv.csv');    % optional extension (public)
files.ff       = fullfile(folders.raw, 'monthly_fed_funds_effective_rate.csv');
files.mort     = fullfile(folders.raw, 'MORTGAGE30US.csv');         % weekly
files.cp       = fullfile(folders.raw, 'RIFSPPFAAD90NB.csv');       % daily
files.tb3      = fullfile(folders.raw, 'DTB3.csv');                 % daily
files.mps_xlsx = fullfile(folders.raw, 'monetary-policy-surprises-data.xlsx');

% Authors (ONLY for FF4 instrument)
files.fac_auth = fullfile(folders.auth15, 'factor_data.csv');

% Output
files.var_out  = fullfile(folders.output, 'VAR_data_extended.csv');
files.fac_out  = fullfile(folders.output, 'Factors_data_extended.csv');

% Key start for GK VAR
startDate = datetime(1979,7,1);

% EBP splice policy:
%  - "none": concatenate GZ then extension with NO level adjustment (default)
%  - "levelshift": level-shift extension to match last overlap mean (OFF by default)
EBP_SPLICE_POLICY = "none";

%% ========================================================================
%  1) LOAD RAW MACRO SERIES (FRED monthly)
% =========================================================================
fprintf('1) Loading raw macro series...\n');

CPI = read_fred_monthly(files.cpi);    % columns: date, value
IP  = read_fred_monthly(files.ip);
GS1 = read_fred_monthly(files.gs1);

CPI.logcpi = 100*log(CPI.value);
IP.logip   = 100*log(IP.value);
GS1.gs1    = GS1.value;

CPI = CPI(:,{'date','logcpi'});
IP  = IP(:, {'date','logip'});
GS1 = GS1(:,{'date','gs1'});

%% ========================================================================
%  2) LOAD YIELDS (GSW) AND BUILD FORWARD (5x5)
% =========================================================================
fprintf('2) Loading GSW yield curve and aggregating to monthly (if needed)...\n');

GSW_daily = read_gsw(files.gsw);              % returns table with date + yields
GSW_m     = gsw_to_monthly_mean(GSW_daily);   % date is first-of-month

% Expected output columns: date, gs2, cm5yr, cm10yr
GSW_m.cm5f5 = (10*GSW_m.cm10yr - 5*GSW_m.cm5yr)/5;

%% ========================================================================
%  3) LOAD EBP (GZ + OPTIONAL EXTENSION)
% =========================================================================
fprintf('3) Loading EBP...\n');

GZ = readtable(files.gz_m);
GZ.date = parse_gz_date(GZ.date); % handles "JAN1973" style
GZ.ebp  = to_num(GZ.ebp);
GZ = GZ(~isnat(GZ.date) & ~isnan(GZ.ebp), :);
GZ = GZ(:,{'date','ebp'});

EBP = GZ;

if exist(files.ebp_ext,'file') == 2
    EXT = readtable(files.ebp_ext);

    % robust parse
    if ~ismember('date', lower(EXT.Properties.VariableNames))
        % assume first column is date
        EXT.date = EXT{:,1};
    end
    if ~isdatetime(EXT.date)
        % try yyyy-MM-dd first, then fallback
        try
            EXT.date = datetime(EXT.date, 'InputFormat','yyyy-MM-dd');
        catch
            EXT.date = datetime(EXT.date);
        end
    end

    % detect ebp column name
    if ~ismember('ebp', lower(EXT.Properties.VariableNames))
        % try second column
        EXT.ebp = EXT{:,2};
    end
    EXT.ebp = to_num(EXT.ebp);
    EXT = EXT(~isnat(EXT.date) & ~isnan(EXT.ebp), :);
    EXT = EXT(:,{'date','ebp'});

    % Concatenate with policy
    EBP = ebp_concatenate(GZ, EXT, EBP_SPLICE_POLICY);
else
    fprintf('   Note: ebp_csv.csv not found; EBP will end at %s\n', string(max(EBP.date)));
end

%% ========================================================================
%  4) LOAD ADDITIONAL SERIES FOR SPREADS (FF, mortgage, CP, TB3)
% =========================================================================
fprintf('4) Loading additional rates for spreads...\n');

FF = read_fred_monthly(files.ff);
FF.ff = FF.value;
FF = FF(:,{'date','ff'});

% Mortgage weekly to monthly mean
M = readtable(files.mort);
M.date = datetime(M.observation_date, 'InputFormat','yyyy-MM-dd');
M.mort30 = to_num(M.MORTGAGE30US);
M = M(~isnan(M.mort30),:);
M.mdate = datetime(year(M.date), month(M.date), 1);
M = groupsummary(M,'mdate','mean','mort30');
M.Properties.VariableNames{'mdate'} = 'date';
M.Properties.VariableNames{'mean_mort30'} = 'mort30';

% CP daily to monthly mean
CP = readtable(files.cp);
CP.date = datetime(CP.observation_date, 'InputFormat','yyyy-MM-dd');
CP.cp90 = to_num(CP.RIFSPPFAAD90NB);
CP = CP(~isnan(CP.cp90),:);
CP.mdate = datetime(year(CP.date), month(CP.date), 1);
CP = groupsummary(CP,'mdate','mean','cp90');
CP.Properties.VariableNames{'mdate'} = 'date';
CP.Properties.VariableNames{'mean_cp90'} = 'cp90';

% TB3 daily -> monthly mean
TB3 = readtable(files.tb3);
TB3.date = datetime(TB3.observation_date, 'InputFormat','yyyy-MM-dd');
TB3.tb3 = to_num(TB3.DTB3);
TB3 = TB3(~isnan(TB3.tb3),:);
TB3.mdate = datetime(year(TB3.date), month(TB3.date), 1);
TB3 = groupsummary(TB3,'mdate','mean','tb3');
TB3.Properties.VariableNames{'mdate'} = 'date';
TB3.Properties.VariableNames{'mean_tb3'} = 'tb3';

%% ========================================================================
%  5) MERGE VAR DATASET
% =========================================================================
fprintf('5) Merging VAR dataset...\n');

T = outerjoin(CPI, IP,  'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, GS1,   'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, GSW_m, 'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, EBP,   'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, FF,    'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, M,     'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, CP,    'Keys','date','MergeKeys',true,'Type','full');
T = outerjoin(T, TB3,   'Keys','date','MergeKeys',true,'Type','full');

% Keep from GK start
T = T(T.date >= startDate, :);

% Construct spreads (pp differences)
T.cp3m_spread_m    = T.cp90 - T.tb3;          % CP - 3m Tbill
T.mortg_spread_m   = T.mort30 - T.cm10yr;     % Mortgage - 10y Treasury

% Create year/month identifiers
T.year  = year(T.date);
T.month = month(T.date);

% Keep columns consistent with authors naming
VAR_out = T(:, {'year','month','logip','logcpi','gs1','ff','gs2','cm5yr','cm10yr','cm5f5','ebp','mortg_spread_m','cp3m_spread_m'});

% Drop rows with missing core variables used in the GK baseline + spreads
core = {'logip','logcpi','gs1','ebp'};
ok = true(height(VAR_out),1);
for k=1:numel(core)
    ok = ok & ~isnan(VAR_out.(core{k}));
end
VAR_out = VAR_out(ok,:);

% Sort
VAR_out = sortrows(VAR_out, {'year','month'});

writetable(VAR_out, files.var_out);
fprintf('   Saved VAR_data_extended.csv with %d months (%s to %s)\n', ...
    height(VAR_out), string(datetime(VAR_out.year(1),VAR_out.month(1),1)), ...
    string(datetime(VAR_out.year(end),VAR_out.month(end),1)));

%% ========================================================================
%  6) BUILD FACTORS (ED2/ED3/ED4 + SP500 monthly using GK Footnote 11 weights)
%      + merge authors FF4 ONLY
%      + build JK(2020) split instruments (PM vs CBI/info)
% =========================================================================
fprintf('6) Building factors...\n');

% Calendar aligned to VAR_out
cal = table(VAR_out.year, VAR_out.month, 'VariableNames', {'year','month'});
cal.date = datetime(cal.year, cal.month, 1);

% Load event-level surprises (Bauer-Swanson file)
FOMC = readtable(files.mps_xlsx, ...
    'Sheet','FOMC (update 2023)', ...
    'VariableNamingRule','preserve');

% ---- Robust detection of date and ED/SP500 columns ----
v = FOMC.Properties.VariableNames;

idxDate = find(contains(lower(v),'date'), 1);
assert(~isempty(idxDate), 'No date column found in FOMC sheet.');
FOMC_date = datetime(FOMC.(v{idxDate}));

idxED2 = find(strcmpi(v,'ED2'),1);
idxED3 = find(strcmpi(v,'ED3'),1);
idxED4 = find(strcmpi(v,'ED4'),1);

assert(~isempty(idxED2), 'No ED2 column found.');
assert(~isempty(idxED3), 'No ED3 column found.');
assert(~isempty(idxED4), 'No ED4 column found.');

idxSP = find(strcmpi(v,'SP500'),1);
if isempty(idxSP)
    idxSP = find(contains(lower(v),'sp500'),1);
end
assert(~isempty(idxSP), 'No SP500 column found (expected SP500 or similar).');

% Footnote-11 style weighting: w=(D-d+1)/D
ED2m   = monthly_weighted_sum(FOMC_date, to_num(FOMC.(v{idxED2})));
ED3m   = monthly_weighted_sum(FOMC_date, to_num(FOMC.(v{idxED3})));
ED4m   = monthly_weighted_sum(FOMC_date, to_num(FOMC.(v{idxED4})));
SP500m = monthly_weighted_sum(FOMC_date, to_num(FOMC.(v{idxSP})));

FAC = cal(:,{'date','year','month'});
FAC = outerjoin(FAC, ED2m,   'Keys','date','MergeKeys',true,'Type','left'); FAC = renamevars(FAC,'value','ed2_tc');
FAC = outerjoin(FAC, ED3m,   'Keys','date','MergeKeys',true,'Type','left'); FAC = renamevars(FAC,'value','ed3_tc');
FAC = outerjoin(FAC, ED4m,   'Keys','date','MergeKeys',true,'Type','left'); FAC = renamevars(FAC,'value','ed4_tc');
FAC = outerjoin(FAC, SP500m, 'Keys','date','MergeKeys',true,'Type','left'); FAC = renamevars(FAC,'value','sp500_tc');

% Months with no meetings -> 0, but only WITHIN each series availability window.
FAC.ed2_tc   = fill_zero_within_window(FAC.date, FAC.ed2_tc,   ED2m.date);
FAC.ed3_tc   = fill_zero_within_window(FAC.date, FAC.ed3_tc,   ED3m.date);
FAC.ed4_tc   = fill_zero_within_window(FAC.date, FAC.ed4_tc,   ED4m.date);
FAC.sp500_tc = fill_zero_within_window(FAC.date, FAC.sp500_tc, SP500m.date);

% Merge authors FF4 only
Afac = readtable(files.fac_auth);
Afac.date = datetime(Afac.year, Afac.month, 1);
Afac = Afac(:,{'date','ff4_tc'});

FAC = outerjoin(FAC, Afac, 'Keys','date','MergeKeys',true,'Type','left');

% Enforce FF4 availability window (paper uses 1991:01–2012:06)
ff4_start = datetime(1991,1,1);
ff4_end   = datetime(2012,6,1);
if ismember('ff4_tc', FAC.Properties.VariableNames)
    FAC.ff4_tc(FAC.date < ff4_start | FAC.date > ff4_end) = NaN;
end

% =========================
% JK(2020) style split instruments
% =========================
sgn_ed2_sp = sign(FAC.ed2_tc) .* sign(FAC.sp500_tc);
FAC.pmneg_ed2_sp500_tc = FAC.ed2_tc .* (sgn_ed2_sp < 0);  % MP-type
FAC.pmpos_ed2_sp500_tc = FAC.ed2_tc .* (sgn_ed2_sp > 0);  % CBI/info-type

% Optional: FF4-based split on the overlap window (will be NaN outside 1991–2012 by construction)
sgn_ff4_sp = sign(FAC.ff4_tc) .* sign(FAC.sp500_tc);
FAC.pmneg_ff4_sp500_tc = FAC.ff4_tc .* (sgn_ff4_sp < 0);
FAC.pmpos_ff4_sp500_tc = FAC.ff4_tc .* (sgn_ff4_sp > 0);

% Output (include the new columns)
FAC_out = FAC(:, {'year','month', ...
    'ff4_tc','ed2_tc','ed3_tc','ed4_tc','sp500_tc', ...
    'pmneg_ed2_sp500_tc','pmpos_ed2_sp500_tc', ...
    'pmneg_ff4_sp500_tc','pmpos_ff4_sp500_tc'});

FAC_out = sortrows(FAC_out, {'year','month'});

writetable(FAC_out, files.fac_out);
fprintf('   Saved Factors_data_extended.csv with %d months\n', height(FAC_out));

fprintf('\nDONE.\n');

%% ========================================================================
%  LOCAL HELPERS
% =========================================================================
function T = read_fred_monthly(path)
    T = readtable(path);
    dcol = T.Properties.VariableNames{1};
    vcol = T.Properties.VariableNames{2};
    T.date  = datetime(T.(dcol));
    T.value = to_num(T.(vcol));
    T = T(:,{'date','value'});
    T = T(~isnat(T.date),:);
    T.date = datetime(year(T.date), month(T.date), 1);
    T = sortrows(T,'date');
    [~,ia] = unique(T.date,'last');
    T = T(ia,:);
end

function X = to_num(x)
    if isnumeric(x); X = x; return; end
    X = str2double(string(x));
end

function G = read_gsw(path)
    G = readtable(path);
    vn = lower(G.Properties.VariableNames);
    idx = find(contains(vn,'date') | contains(vn,'observation_date'),1);
    assert(~isempty(idx), 'GSW file: could not find a date column');
    dname = G.Properties.VariableNames{idx};
    G.date = datetime(G.(dname));
    G = G(~isnat(G.date),:);

    if ismember('SVENY02', G.Properties.VariableNames)
        y2  = to_num(G.SVENY02);
        y5  = to_num(G.SVENY05);
        y10 = to_num(G.SVENY10);
    elseif any(strcmpi(G.Properties.VariableNames,'gs2'))
        y2  = to_num(G{:, find(strcmpi(G.Properties.VariableNames,'gs2'))});
        y5  = to_num(G{:, find(strcmpi(G.Properties.VariableNames,'cm5yr'))});
        y10 = to_num(G{:, find(strcmpi(G.Properties.VariableNames,'cm10yr'))});
    else
        error('GSW file: unsupported column naming. Need SVENY02/SVENY05/SVENY10 or gs2/cm5yr/cm10yr');
    end

    G = table(G.date, y2, y5, y10, 'VariableNames', {'date','gs2','cm5yr','cm10yr'});
end

function M = gsw_to_monthly_mean(G)
    G.mdate = datetime(year(G.date), month(G.date), 1);
    M = groupsummary(G,'mdate','mean',{'gs2','cm5yr','cm10yr'});
    M.Properties.VariableNames{'mdate'} = 'date';
    M.Properties.VariableNames{'mean_gs2'} = 'gs2';
    M.Properties.VariableNames{'mean_cm5yr'} = 'cm5yr';
    M.Properties.VariableNames{'mean_cm10yr'} = 'cm10yr';
end

function out = monthly_weighted_sum(dates, shock)
    ok = ~isnat(dates) & ~isnan(shock);
    dates = dates(ok); shock = shock(ok);

    mdate = datetime(year(dates), month(dates), 1);
    d = day(dates);
    D = eomday(year(dates), month(dates));
    w = (D - d + 1) ./ D;

    contrib = shock .* w;

    [g, mvals] = findgroups(mdate);
    z = splitapply(@sum, contrib, g);

    out = table(mvals, z, 'VariableNames', {'date','value'});
end

function x = fill_zero_within_window(all_dates, x, avail_dates)
    % Set NaN -> 0 only within [min(avail_dates), max(avail_dates)].
    % Outside that window, keep NaN (series not available).
    if isempty(avail_dates)
        return;
    end
    w0 = min(avail_dates);
    w1 = max(avail_dates);
    in = (all_dates >= w0) & (all_dates <= w1);
    x(in & isnan(x)) = 0;
end

function d = parse_gz_date(x)
    s = string(x);
    s = strtrim(s);
    mon = extractBetween(upper(s), 1, 3);
    yr  = extractAfter(s, 3);
    d = datetime(strcat("01-", mon, "-", yr), 'InputFormat','dd-MMM-yyyy', 'Locale','en_US');
end

function E = ebp_concatenate(GZ, EXT, policy)
    GZ  = sortrows(GZ,'date');
    EXT = sortrows(EXT,'date');

    GZ.date  = datetime(year(GZ.date),  month(GZ.date),  1);
    EXT.date = datetime(year(EXT.date), month(EXT.date), 1);

    overlap_start = max(min(GZ.date), min(EXT.date));
    overlap_end   = min(max(GZ.date), max(EXT.date));

    if overlap_start <= overlap_end
        a = GZ(GZ.date>=overlap_start & GZ.date<=overlap_end, :);
        b = EXT(EXT.date>=overlap_start & EXT.date<=overlap_end, :);

        a.ebp_gz = a.ebp;  a.ebp = [];
        b.ebp_ext = b.ebp; b.ebp = [];

        J = innerjoin(a, b, 'Keys', 'date');  % no Suffixes (compat)

        if ~isempty(J)
            gap = mean(J.ebp_gz - J.ebp_ext, 'omitnan');
            fprintf('   EBP overlap %s–%s: mean(GZ-EXT)=%.4f pp\n', ...
                string(overlap_start), string(overlap_end), gap);

            if policy == "levelshift"
                EXT.ebp = EXT.ebp + gap;
                fprintf('   Applied LEVEL SHIFT to EXT by %.4f pp (policy=levelshift)\n', gap);
            else
                fprintf('   No level adjustment applied (policy=none)\n');
            end
        end
    else
        fprintf('   No overlap between GZ and EXT EBP; concatenating as-is.\n');
    end

    cutoff = min(EXT.date);
    E = [GZ(GZ.date < cutoff,:); EXT];
end

