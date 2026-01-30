% VAR_main_Extension.m
% =========================================================================

clearvars; close all; clc;

fprintf('=========================================================================\n');
fprintf('STEP 2: REPLICATING GERTLER & KARADI (2015) FIGURE 1\n');
fprintf('=========================================================================\n\n');

%% Import and load extended time-series data
fprintf('Step 1: Loading data...\n');
Import_Data;  % Must create DATASET_VAR and DATASET_FACTORS in base workspace
fprintf('  Data loaded from extended datasets\n\n');

%% Bootstrap and display settings
nboot  = 1000;   % Paper uses 10000; 1000 is fine for development
clevel = 95;     % Percentile bands
VAR.fontsize = 14;

VAR.switch_extern = 0;   % Keep simple for Figure 1
VAR.switch_exp    = 1;   % In original code used for expectations VARs

%% ========================================================================
% STEP 2: REPLICATION SETTINGS (1979-2012)
% =========================================================================
fprintf('Step 2: Setting replication parameters...\n');

% Sample period (frozen for Step 2 replication)
smpl_min_VAR_vec = [1979 7];
smpl_max_VAR_vec = [2023 12];

fprintf('  Sample period: %d:%02d to %d:%02d (REPLICATION)\n', ...
    smpl_min_VAR_vec(1), smpl_min_VAR_vec(2), smpl_max_VAR_vec(1), smpl_max_VAR_vec(2));

% Monetary policy variables
monpol_vars_cell              = {'GS1'};
monpol_vars_label_cell        = {'1 year rate'};
monpol_vars_label_short_cell  = {'1YR'};

% Factors (instruments) by policy indicator
% Fed Funds policy indicator case (not used in baseline Figure 1 here)
factors1_cell_FF       = {'MP1_TC'};
factors1_label_cell_FF = {'MP1_TC'};
factors2_cell_FF       = {''};

% GS1 policy indicator case (baseline in GK Figure 1 uses FF4)
factors1_cell_GS1       = {'ED2_TC'};
factors1_label_cell_GS1 = {'ED2_TC'};
factors2_cell_GS1       = {''};

% GS2 policy indicator case (not used here, but kept for completeness)
factors1_cell_GS2       = {'MP1_TC'};
factors1_label_cell_GS2 = {'GSS'};
factors2_cell_GS2       = {'FF4_TC','ED2_TC','ED3_TC','ED4_TC'};

% Factor availability window (replication)
smpl_min_factors_vec = [1991 1];
smpl_max_factors_vec = [2023 12];

fprintf('  Instrument window (intended): %d:%02d to %d:%02d\n', ...
    smpl_min_factors_vec(1), smpl_min_factors_vec(2), smpl_max_factors_vec(1), smpl_max_factors_vec(2));

figure_name = 'EXTENSION';

%% Credit spreads (for Figure 1)
spreads_cell             = {'EBP','MORTG_SPREAD','CP3M_SPREAD'};
spreads_label_cell       = {'Excess Bond Premium','Mortgage spread','Commercial Paper spread (3 months)'};
spreads_label_short_cell = {'EBP','MORTG.','CP3M'};
no_spread_single         = 1;   % Figure 1 uses EBP

%% External variables (kept for compatibility; not used in Figure 1 code path)
extern_vars_cell              = {'FF','GS1','GS2','CM5YR','CM10YR','CM5F5','FF_EXP1YR'};
extern_vars_label_cell        = {'Federal Funds Rate','1 year rate','2 year rate','5 year rate','10 year rate','5x5 forward','1 year expectations (FF)'};
extern_vars_label_short_cell  = {'FF','1YR','2YR','5YR','10YR','5F5','EXP1YR(FF)'};
extern_vars_smpl_min_vec      = [ones(6,1)*[1979 7]; [1983 3]];

%% VAR settings
VAR.irhor = 48;
VAR.p     = 12;

fprintf('  VAR lags: %d\n', VAR.p);
fprintf('  IRF horizon: %d months\n\n', VAR.irhor);

%% Count variables
no_smpl_max_VAR  = size(smpl_max_VAR_vec,1);
no_monpol_vars   = numel(monpol_vars_cell);
no_spreads       = numel(no_spread_single);

fprintf('Step 3: Data structures ready (from Import_Data.m)\n\n');

%% ========================================================================
% RUN VAR ESTIMATION
% =========================================================================
fprintf('Step 4: Running VAR estimation and IRF computation...\n');
fprintf('  This may take several minutes depending on bootstrap replications...\n\n');

ii = cell(no_monpol_vars,1);

for ii_monpol = 1:no_monpol_vars

    % Fed funds indicator flag (kept for GK code logic)
    if strcmp(monpol_vars_cell{ii_monpol},'FF')
        VAR.monpol_FF = 'yes';
    else
        VAR.monpol_FF = 'no';
    end

    ii{ii_monpol} = 1000*(ii_monpol-1)+1;

    % Choose factor sets based on policy indicator
    if strcmp(monpol_vars_cell{ii_monpol},'GS1')
        factors1_cell       = factors1_cell_GS1;
        factors1_label_cell = factors1_label_cell_GS1;
        factors2_cell       = factors2_cell_GS1;
    elseif strcmp(monpol_vars_cell{ii_monpol},'FF')
        factors1_cell       = factors1_cell_FF;
        factors1_label_cell = factors1_label_cell_FF;
        factors2_cell       = factors2_cell_FF;
    elseif strcmp(monpol_vars_cell{ii_monpol},'GS2')
        factors1_cell       = factors1_cell_GS2;
        factors1_label_cell = factors1_label_cell_GS2;
        factors2_cell       = factors2_cell_GS2;
    else
        error('Unknown monpol var: %s', monpol_vars_cell{ii_monpol});
    end

    no_factors1 = numel(factors1_cell);
    no_factors2 = numel(factors2_cell);

    % For replication we typically have just one smpl_max_VAR_vec row
    for ii_smpl_max = 1:no_smpl_max_VAR
        for ii_factors = 1:no_factors1
            for ii_spreads = 1:no_spreads

                % ----------------------------
                % Sample windows (safe datetime arithmetic)
                % ----------------------------
                smpl_min_VAR = smpl_min_VAR_vec(1,:);
                smpl_max_VAR = smpl_max_VAR_vec(ii_smpl_max,:);

                varStartDT = datetime(smpl_min_VAR(1,1), smpl_min_VAR(1,2), 1);
                varEndDT   = datetime(smpl_max_VAR(1,1), smpl_max_VAR(1,2), 1);

                instStartDT = datetime(smpl_min_factors_vec(1,1), smpl_min_factors_vec(1,2), 1);
                instEndDT   = datetime(smpl_max_factors_vec(1,1), smpl_max_factors_vec(1,2), 1);

                % Factors start at max(instrument start, VAR start + p months)
                lagStartDT  = varStartDT + calmonths(VAR.p);
                smplMinFDT  = max(instStartDT, lagStartDT);

                % Factors end at min(VAR end, instrument end)
                smplMaxFDT  = min(varEndDT, instEndDT);

                smpl_min_FACTORS = [year(smplMinFDT), month(smplMinFDT)];
                smpl_max_FACTORS = [year(smplMaxFDT), month(smplMaxFDT)];

                % ----------------------------
                % Logging
                % ----------------------------
                if strcmp(factors2_cell{1,1},'')
                    fac2label = '';
                else
                    fac2label = [' + ' strjoin(factors2_cell, ',')];
                end

                fprintf('\n#%3.0f\n', ii{ii_monpol});
                fprintf('MONPOL:  %s\n', monpol_vars_cell{ii_monpol});
                fprintf('SPREAD:  %s\n', spreads_cell{no_spread_single(ii_spreads)});
                fprintf('FACTOR:  %s%s\n', factors1_label_cell{ii_factors}, fac2label);
                fprintf('SAMPLE:  %d-%02d to %d-%02d\n', smpl_min_VAR(1,1), smpl_min_VAR(1,2), smpl_max_VAR(1,1), smpl_max_VAR(1,2));
                fprintf('FACTORS: %d-%02d to %d-%02d\n', smpl_min_FACTORS(1,1), smpl_min_FACTORS(1,2), smpl_max_FACTORS(1,1), smpl_max_FACTORS(1,2));

                % ----------------------------
                % Date indexing with uniqueness assertions
                % ----------------------------
                VAR.smpl_min_VAR = local_unique_index(DATASET_VAR, smpl_min_VAR(1,1), smpl_min_VAR(1,2), 'VAR start');
                VAR.smpl_max_VAR = local_unique_index(DATASET_VAR, smpl_max_VAR(1,1), smpl_max_VAR(1,2), 'VAR end');

                VAR.smpl_max_VAR_factors = local_unique_index(DATASET_VAR, smpl_max_FACTORS(1,1), smpl_max_FACTORS(1,2), 'VAR factors-end anchor');

                VAR.smpl_min_FACTORS = local_unique_index(DATASET_FACTORS, smpl_min_FACTORS(1,1), smpl_min_FACTORS(1,2), 'FACTORS start');
                VAR.smpl_max_FACTORS = local_unique_index(DATASET_FACTORS, smpl_max_FACTORS(1,1), smpl_max_FACTORS(1,2), 'FACTORS end');

                % ----------------------------
                % Variable selection (Proxy-SVAR: MonPol must be first)
                % Ordering: [MonPol, CPI, IP, Spread]
                % ----------------------------
                VAR.select_vars = {monpol_vars_cell{ii_monpol}};
                VAR.select_vars_label = {monpol_vars_label_cell{ii_monpol}};
                VAR.select_vars_label_short = {monpol_vars_label_short_cell{ii_monpol}};

                VAR.select_vars = [VAR.select_vars, {'LCPI','LIP'}];
                VAR.select_vars_label = [VAR.select_vars_label, {'CPI','IP'}];
                VAR.select_vars_label_short = [VAR.select_vars_label_short, {'CPI','IP'}];

                spreadName = spreads_cell{no_spread_single(ii_spreads)};
                VAR.select_vars = [VAR.select_vars, {spreadName}];
                VAR.select_vars_label = [VAR.select_vars_label, {spreads_label_cell{no_spread_single(ii_spreads)}}];
                VAR.select_vars_label_short = [VAR.select_vars_label_short, {spreads_label_short_cell{no_spread_single(ii_spreads)}}];

                % Cholesky VAR order: [CPI, IP, MonPol, Spread]
                VAR.chol_order = [2,3,1,4];

                % ----------------------------
                % Factor selection (labels = names; robust when multiple instruments)
                % ----------------------------
                VAR.select_factors = {factors1_cell{ii_factors}};
                VAR.select_factors_label = {factors1_label_cell{ii_factors}};

                if ~strcmp(factors2_cell{1,1},'')
                    for jj_factors2 = 1:no_factors2
                        VAR.select_factors = [VAR.select_factors, {factors2_cell{jj_factors2}}];
                        VAR.select_factors_label = [VAR.select_factors_label, {factors2_cell{jj_factors2}}];
                    end
                end

                % ----------------------------
                % Data extraction
                % ----------------------------
                no_vars_VAR = numel(VAR.select_vars);
                VAR.vars = zeros(VAR.smpl_max_VAR, no_vars_VAR);
                for jj = 1:no_vars_VAR
                    VAR.vars(:,jj) = DATASET_VAR.TSERIES(1:VAR.smpl_max_VAR, DATASET_VAR.MAP(VAR.select_vars{jj}));
                end
                VAR.vars = VAR.vars(VAR.smpl_min_VAR:VAR.smpl_max_VAR,:);

                no_factors_selected = numel(VAR.select_factors);
                Tfac = VAR.smpl_max_FACTORS - VAR.smpl_min_FACTORS + 1;
                VAR.proxies = zeros(Tfac, no_factors_selected);
                for jj = 1:no_factors_selected
                    VAR.proxies(:,jj) = DATASET_FACTORS.TSERIES(VAR.smpl_min_FACTORS:VAR.smpl_max_FACTORS, ...
                        DATASET_FACTORS.MAP(VAR.select_factors{jj}));
                end

                % Diagnostics
                fprintf('  Proxy extraction:\n');
                fprintf('    Factor(1): %s\n', VAR.select_factors{1});
                fprintf('    Obs: %d\n', size(VAR.proxies,1));
                fprintf('    Non-NaN: %d\n', sum(~isnan(VAR.proxies(:,1))));
                fprintf('    Mean: %.4f | Std: %.4f\n', nanmean(VAR.proxies(:,1)), nanstd(VAR.proxies(:,1)));

                % Extra guard: FF4 should not be artificially zero-filled
                if strcmpi(VAR.select_factors{1}, 'FF4_TC')
                    shareZero = mean(VAR.proxies(:,1)==0, 'omitnan');
                    fprintf('    FF4 zero-share (within factor sample) = %.3f\n', shareZero);
                end

                % ----------------------------
                % Defaults to keep GK code stable for Figure 1
                % ----------------------------
                VAR.T_m_end = 0;
                VAR.term_spreads = [];
                VAR.term_spreads_matur = 0;
                VAR.term_spreads_init = [];
                VAR.real_rates = [];
                VAR.real_rates_init = 0;
                VAR.real_rates_matur = 0;
                VAR.excess_return = [];
                VAR.excess_return_matur = [];

                % ----------------------------
                % Output path (relative to this file)
                % ----------------------------
                rootDir = fileparts(mfilename('fullpath'));
                output_dir = fullfile(rootDir, 'output');
                if ~exist(output_dir, 'dir')
                    mkdir(output_dir);
                    fprintf('  Created output directory: %s\n', output_dir);
                end
                VAR.figure_name = fullfile(output_dir, figure_name);

                % ----------------------------
                % Proxy SVAR + bootstrap
                % ----------------------------
                fprintf('  Estimating Proxy SVAR...\n');
                VAR = doProxySVAR_single(VAR);
                fprintf('  IRFs computed successfully\n');

                fprintf('  Running bootstrap (%d replications)...\n', nboot);
                VARbs = doProxySVARbootstrap_single(VAR, nboot, clevel);

                VAR.irsH = VARbs.irsH;
                VAR.irsL = VARbs.irsL;
                if isfield(VARbs, 'irsH_e')
                    VAR.irsH_e = VARbs.irsH_e;
                    VAR.irsL_e = VARbs.irsL_e;
                end

                % ----------------------------
                % Cholesky comparison + bootstrap
                % ----------------------------
                fprintf('  Estimating Cholesky SVAR for comparison...\n');

                VAR.k = 1;
                VARChol = VAR;
                VARChol.vars = VAR.vars(:, VAR.chol_order);

                VARChol = doCholSVAR_single(VARChol);
                VARCholbs = doCholSVARbootstrap_single(VARChol, nboot, clevel);

                % ----------------------------
                % Plot (side-by-side)
                % ----------------------------
                fprintf('  Generating figures...\n');

                figure_num = ii{ii_monpol} + ii_smpl_max*10 + ii_factors*100 + ii_spreads;
                nRow = VAR.n;   % 4 rows
                nCol = 2;       % EI left / Chol right

                plot_figure_sep(VAR, VARChol, VARbs, VARCholbs, nCol, nRow, figure_num, VAR.switch_extern);
                fprintf('  Figure generation complete\n\n');

            end
        end
    end
end

fprintf('=========================================================================\n');
fprintf('STEP 2 COMPLETE: Figure 1 should now be saved in the output/ folder\n');
fprintf('=========================================================================\n\n');

fprintf('For Step 3 (Extension):\n');
fprintf('  Change smpl_max_VAR_vec to [2023 12] (or latest)\n');
fprintf('  Change smpl_max_factors_vec to [2023 12] (or latest)\n');
fprintf('  Switch factors1_cell_GS1 to {''ED2_TC''} (or chosen proxy)\n');
fprintf('  Set figure_name = ''EXTENSION''\n\n');

fprintf('For Step 4 (Extension Pre-Covid):\n');
fprintf('  Change smpl_max_VAR_vec to [2019 12]\n');
fprintf('  Change smpl_max_factors_vec to [2019 12]\n');
fprintf('  Switch factors1_cell_GS1 to {''ED2_TC''}\n');
fprintf('  Set figure_name = ''EXTENSION_PRECOVID''\n\n');

%% ========================================================================
% Local helper: unique date index in DATASET struct
% =========================================================================
function idx = local_unique_index(DATASET, yy, mm, label)
    YEAR  = DATASET.TSERIES(:, DATASET.MAP('YEAR'));
    MONTH = DATASET.TSERIES(:, DATASET.MAP('MONTH'));
    idx = find(YEAR==yy & MONTH==mm);
    assert(numel(idx)==1, '%s date %d-%02d appears %d times (expected exactly 1).', label, yy, mm, numel(idx));
end