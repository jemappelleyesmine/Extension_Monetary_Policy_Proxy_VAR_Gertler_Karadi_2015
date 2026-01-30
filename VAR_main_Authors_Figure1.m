% VAR_main_Authors_Figure1.m
% =========================================================================

clearvars; close all; clc;

fprintf('=========================================================================\n');
fprintf('AUTHORS DATA RUN: GERTLER & KARADI (2015) FIGURE 1\n');
fprintf('=========================================================================\n\n');

%% Load AUTHORS datasets
fprintf('Step 1: Loading AUTHORS data...\n');
Import_Data_Authors;
fprintf('  Authors data loaded.\n\n');

%% Bootstrap and display settings
nboot  = 1000;   % Paper uses 10000; set higher later if needed
clevel = 95;
VAR.fontsize = 14;

VAR.switch_extern = 0;   % Figure 1
VAR.switch_exp    = 1;

%% Replication sample (as in paper)
fprintf('Step 2: Setting replication parameters...\n');

smpl_min_VAR_vec = [1979 7];
smpl_max_VAR_vec = [2012 6];

fprintf('  Sample period: %d:%02d to %d:%02d\n', ...
    smpl_min_VAR_vec(1), smpl_min_VAR_vec(2), smpl_max_VAR_vec(1), smpl_max_VAR_vec(2));

% Policy indicator
monpol_vars_cell             = {'GS1'};
monpol_vars_label_cell       = {'1 year rate'};
monpol_vars_label_short_cell = {'1YR'};

% Factors: baseline uses FF4 for GS1
factors1_cell_GS1       = {'FF4_TC'};
factors1_label_cell_GS1 = {'FF4_TC'};
factors2_cell_GS1       = {''};

% Instrument availability window
smpl_min_factors_vec = [1991 1];
smpl_max_factors_vec = [2012 6];

fprintf('  Instrument window: %d:%02d to %d:%02d\n', ...
    smpl_min_factors_vec(1), smpl_min_factors_vec(2), smpl_max_factors_vec(1), smpl_max_factors_vec(2));

% Output name
figure_name = 'AUTHORS_FIGURE1';

%% Figure 1 spread: EBP
spreads_cell             = {'EBP','MORTG_SPREAD','CP3M_SPREAD'};
spreads_label_cell       = {'Excess Bond Premium','Mortgage spread','Commercial Paper spread (3 months)'};
spreads_label_short_cell = {'EBP','MORTG.','CP3M'};
no_spread_single         = 1;  % EBP

%% VAR settings
VAR.irhor = 48;
VAR.p     = 12;

fprintf('  VAR lags: %d\n', VAR.p);
fprintf('  IRF horizon: %d months\n\n', VAR.irhor);

%% RUN
fprintf('Step 3: Running VAR + Proxy SVAR + bootstrap...\n');

% Choose factors for GS1 case
factors1_cell       = factors1_cell_GS1;
factors1_label_cell = factors1_label_cell_GS1;
factors2_cell       = factors2_cell_GS1;

% Safe datetime sample arithmetic
smpl_min_VAR = smpl_min_VAR_vec;
smpl_max_VAR = smpl_max_VAR_vec;

varStartDT = datetime(smpl_min_VAR(1), smpl_min_VAR(2), 1);
varEndDT   = datetime(smpl_max_VAR(1), smpl_max_VAR(2), 1);

instStartDT = datetime(smpl_min_factors_vec(1), smpl_min_factors_vec(2), 1);
instEndDT   = datetime(smpl_max_factors_vec(1), smpl_max_factors_vec(2), 1);

lagStartDT  = varStartDT + calmonths(VAR.p);
smplMinFDT  = max(instStartDT, lagStartDT);
smplMaxFDT  = min(varEndDT, instEndDT);

smpl_min_FACTORS = [year(smplMinFDT), month(smplMinFDT)];
smpl_max_FACTORS = [year(smplMaxFDT), month(smplMaxFDT)];

fprintf('  VAR sample:     %d-%02d to %d-%02d\n', smpl_min_VAR(1), smpl_min_VAR(2), smpl_max_VAR(1), smpl_max_VAR(2));
fprintf('  Factors sample: %d-%02d to %d-%02d\n\n', smpl_min_FACTORS(1), smpl_min_FACTORS(2), smpl_max_FACTORS(1), smpl_max_FACTORS(2));

% Unique date indices
VAR.smpl_min_VAR = local_unique_index(DATASET_VAR, smpl_min_VAR(1), smpl_min_VAR(2), 'VAR start');
VAR.smpl_max_VAR = local_unique_index(DATASET_VAR, smpl_max_VAR(1), smpl_max_VAR(2), 'VAR end');

VAR.smpl_min_FACTORS = local_unique_index(DATASET_FACTORS, smpl_min_FACTORS(1), smpl_min_FACTORS(2), 'FACTORS start');
VAR.smpl_max_FACTORS = local_unique_index(DATASET_FACTORS, smpl_max_FACTORS(1), smpl_max_FACTORS(2), 'FACTORS end');

% Select VAR variables: [MonPol, CPI, IP, Spread]
VAR.select_vars = {monpol_vars_cell{1}};
VAR.select_vars_label = {monpol_vars_label_cell{1}};
VAR.select_vars_label_short = {monpol_vars_label_short_cell{1}};

VAR.select_vars = [VAR.select_vars, {'LCPI','LIP'}];
VAR.select_vars_label = [VAR.select_vars_label, {'CPI','IP'}];
VAR.select_vars_label_short = [VAR.select_vars_label_short, {'CPI','IP'}];

spreadName = spreads_cell{no_spread_single};
VAR.select_vars = [VAR.select_vars, {spreadName}];
VAR.select_vars_label = [VAR.select_vars_label, {spreads_label_cell{no_spread_single}}];
VAR.select_vars_label_short = [VAR.select_vars_label_short, {spreads_label_short_cell{no_spread_single}}];

% Cholesky order: [CPI, IP, MonPol, Spread]
VAR.chol_order = [2,3,1,4];

% Factors
VAR.select_factors = {factors1_cell{1}};
VAR.select_factors_label = {factors1_label_cell{1}};
if ~strcmp(factors2_cell{1},'')
    for jj = 1:numel(factors2_cell)
        VAR.select_factors = [VAR.select_factors, {factors2_cell{jj}}];
        VAR.select_factors_label = [VAR.select_factors_label, {factors2_cell{jj}}];
    end
end

% Extract VAR variables
no_vars_VAR = numel(VAR.select_vars);
VAR.vars = zeros(VAR.smpl_max_VAR, no_vars_VAR);
for jj = 1:no_vars_VAR
    VAR.vars(:,jj) = DATASET_VAR.TSERIES(1:VAR.smpl_max_VAR, DATASET_VAR.MAP(VAR.select_vars{jj}));
end
VAR.vars = VAR.vars(VAR.smpl_min_VAR:VAR.smpl_max_VAR,:);

% Extract proxies
no_factors_selected = numel(VAR.select_factors);
Tfac = VAR.smpl_max_FACTORS - VAR.smpl_min_FACTORS + 1;
VAR.proxies = zeros(Tfac, no_factors_selected);
for jj = 1:no_factors_selected
    VAR.proxies(:,jj) = DATASET_FACTORS.TSERIES(VAR.smpl_min_FACTORS:VAR.smpl_max_FACTORS, ...
        DATASET_FACTORS.MAP(VAR.select_factors{jj}));
end

fprintf('  Proxy(1) = %s | Obs=%d | Non-NaN=%d | Mean=%.4f | Std=%.4f\n\n', ...
    VAR.select_factors{1}, size(VAR.proxies,1), sum(~isnan(VAR.proxies(:,1))), ...
    nanmean(VAR.proxies(:,1)), nanstd(VAR.proxies(:,1)));

% Defaults
VAR.T_m_end = 0;
VAR.term_spreads = [];
VAR.term_spreads_matur = 0;
VAR.term_spreads_init = [];
VAR.real_rates = [];
VAR.real_rates_init = 0;
VAR.real_rates_matur = 0;
VAR.excess_return = [];
VAR.excess_return_matur = [];

% Output path relative to this file
rootDir = fileparts(mfilename('fullpath'));
output_dir = fullfile(rootDir, 'output');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
VAR.figure_name = fullfile(output_dir, figure_name);

% Estimate proxy SVAR + bootstrap
VAR = doProxySVAR_single(VAR);
VARbs = doProxySVARbootstrap_single(VAR, nboot, clevel);
VAR.irsH = VARbs.irsH; VAR.irsL = VARbs.irsL;
if isfield(VARbs,'irsH_e'); VAR.irsH_e = VARbs.irsH_e; VAR.irsL_e = VARbs.irsL_e; end

% Cholesky + bootstrap
VAR.k = 1;
VARChol = VAR;
VARChol.vars = VAR.vars(:, VAR.chol_order);

VARChol = doCholSVAR_single(VARChol);
VARCholbs = doCholSVARbootstrap_single(VARChol, nboot, clevel);

% Plot
figure_num = 100; % arbitrary, keeps separate from other runs
nRow = VAR.n; nCol = 2;
plot_figure_sep(VAR, VARChol, VARbs, VARCholbs, nCol, nRow, figure_num, VAR.switch_extern);

fprintf('DONE. Authors Figure 1 output should be in output/ with prefix: %s\n', figure_name);

%% Helper
function idx = local_unique_index(DATASET, yy, mm, label)
    YEAR  = DATASET.TSERIES(:, DATASET.MAP('YEAR'));
    MONTH = DATASET.TSERIES(:, DATASET.MAP('MONTH'));
    idx = find(YEAR==yy & MONTH==mm);
    assert(numel(idx)==1, '%s date %d-%02d appears %d times (expected exactly 1).', label, yy, mm, numel(idx));
end