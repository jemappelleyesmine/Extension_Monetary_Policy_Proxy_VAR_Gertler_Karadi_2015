% =========================================================================
% JK (2020) DECOMPOSITION ON EXTENDED DATA
% =========================================================================

clear all; close all;

%% ========================================================================
% CONFIGURATION
% =========================================================================

% Which analysis to run
RUN_BASELINE = true;   % JK baseline sign restrictions (sgnm2)
RUN_POORMAN = true;    % JK poor man's (cholesky on proxies)

% Sample periods
% NOTE: ED2/SP500 instruments available from 1988-02, so PreCovid/Extension 
%       can start earlier than Replication period for maximum data usage
samples = {
    'Replication', [1991  1; 2012  6];   % JK (2020) replication period
    'PreCovid',    [1988  2; 2019 12];   % Maximum pre-COVID (ED2/SP500 available from 1988-02)
    'Extension',   [1988  2; 2023 12]    % Full extension with COVID
};

% Data paths
path_var = 'data/output_data/VAR_data_extended.csv';
path_fac = 'data/output_data/Factors_data_extended.csv';

% Output directory
output_dir = 'output';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% VAR settings
prior.lags = 12;
prior.minnesota.tightness = 0.2;
prior.minnesota.decay = 1;

% Gibbs sampler settings
gssettings.ndraws = 1000;   % Increase to 4000 for final results
gssettings.burnin = 1000;   % Increase to 4000 for final results  
gssettings.saveevery = 1;
gssettings.computemarglik = 0;

% IRF horizon
MAlags = 48;

%% ========================================================================
% LOAD AND PREPARE DATA
% =========================================================================

fprintf('Loading extended datasets...\n');

% Load CSVs
VAR_tbl = readtable(path_var);
FAC_tbl = readtable(path_fac);

% Validate required columns exist
required_var_cols = {'year', 'month', 'logip', 'logcpi', 'gs1', 'ebp'};
required_fac_cols = {'year', 'month', 'ed2_tc', 'sp500_tc', 'pmneg_ed2_sp500_tc', 'pmpos_ed2_sp500_tc'};

missing_var = setdiff(required_var_cols, VAR_tbl.Properties.VariableNames);
missing_fac = setdiff(required_fac_cols, FAC_tbl.Properties.VariableNames);

if ~isempty(missing_var)
    error('VAR_data_extended.csv missing columns: %s', strjoin(missing_var, ', '));
end
if ~isempty(missing_fac)
    error('Factors_data_extended.csv missing columns: %s', strjoin(missing_fac, ', '));
end

fprintf('  VAR_tbl: %d rows, columns: %s\n', height(VAR_tbl), strjoin(VAR_tbl.Properties.VariableNames, ', '));
fprintf('  FAC_tbl: %d rows, columns: %s\n', height(FAC_tbl), strjoin(FAC_tbl.Properties.VariableNames, ', '));

% Join on year/month directly (avoids duplicated year/month columns)
merged = innerjoin(VAR_tbl, FAC_tbl, 'Keys', {'year','month'});

% Create yyyymm (for clean filtering)
merged.yyyymm = merged.year*100 + merged.month;

fprintf('  Merged data: %d observations (%d-%02d to %d-%02d)\n', ...
    height(merged), merged.year(1), merged.month(1), ...
    merged.year(end), merged.month(end));
fprintf('  Merged columns: %s\n', strjoin(merged.Properties.VariableNames, ', '));

%% ========================================================================
% LOOP OVER SAMPLES AND METHODOLOGIES
% =========================================================================

for s = 1:size(samples, 1)
    sample_name = samples{s, 1};
    spl = samples{s, 2};
    
    fprintf('\n========================================\n');
    fprintf('SAMPLE: %s (%d-%02d to %d-%02d)\n', sample_name, spl(1,1), spl(1,2), spl(2,1), spl(2,2));
    fprintf('========================================\n');
    
    % Filter data to sample period (use yyyymm for correct date comparison)
    start_yyyymm = spl(1,1)*100 + spl(1,2);
    end_yyyymm   = spl(2,1)*100 + spl(2,2);
    
    idx = (merged.yyyymm >= start_yyyymm) & (merged.yyyymm <= end_yyyymm);
    sample_data = merged(idx, :);
    
    %% --------------------------------------------------------------------
    % BASELINE: Sign restrictions on (ED2, SP500)
    % --------------------------------------------------------------------
    if RUN_BASELINE
        fprintf('\n--- BASELINE: Sign restrictions ---\n');
        
        % Prepare data structure (JK format)
        data = struct();
        
        % High-frequency instruments (i.i.d. block)
        data.y = [sample_data.ed2_tc, ...           % m1: ED2 surprise
                  sample_data.sp500_tc, ...         % m2: SP500 surprise
                  sample_data.gs1, ...              % y1: 1-year rate
                  sample_data.logcpi, ...           % y2: log CPI
                  sample_data.logip, ...            % y3: log IP
                  sample_data.ebp];                 % y4: EBP
        
        data.names = {'ed2_tc', 'sp500_tc', 'gs1', 'logcpi', 'logip', 'ebp'};
        data.Nm = 2;  % First 2 variables are i.i.d. instruments
        data.w = ones(size(data.y, 1), 1);  % Constant term
        
        % Time vector (JK format)
        ym2t = @(y, m) y + m/12 - 1/24;
        data.time = ym2t(sample_data.year, sample_data.month);
        
        % Minnesota prior settings
        prior.Nm = data.Nm;
        prior.minnesota.mvector = [0, 0, 1, 0, 0, 1];  % Unit roots for nonstationary vars
        
        % Check data quality
        fprintf('  Data dimensions: T=%d, N=%d, Nm=%d\n', ...
            size(data.y,1), length(data.names), data.Nm);
        fprintf('  ED2  NaNs: %d (%.1f%%)\n', sum(isnan(data.y(:,1))), 100*mean(isnan(data.y(:,1))));
        fprintf('  SP500 NaNs: %d (%.1f%%)\n', sum(isnan(data.y(:,2))), 100*mean(isnan(data.y(:,2))));
        
        % Estimate VAR with i.i.d. instruments
        fprintf('  Estimating VAR... ');
        res_baseline = VAR_withiid1kf(data, prior, gssettings);
        fprintf('Done!\n');
        
        % Apply sign restrictions
        fprintf('  Applying sign restrictions... ');
        N = length(data.names);
        dims = {[1 2]};  % Rotate first two shocks (ED2, SP500)
        imonpol = 1; inews = 2;
        
        % Sign restrictions:
        % MP shock: ED2↑, SP500↓ (contractionary)
        % CBI shock: ED2↑, SP500↑ (information)
        test_restr = @(irfs) ...
            irfs(1,imonpol,1) > 0 && irfs(2,imonpol,1) < 0 && ...  % MP
            irfs(1,inews,1) > 0 && irfs(2,inews,1) > 0;            % CBI
        
        b_normalize = ones(1, N);
        max_try = 1000;
        
        irfs_baseline = resirfssign(res_baseline, MAlags, dims, test_restr, b_normalize, max_try);
        fprintf('Done!\n');
        
        % Store results
        results_baseline.(sample_name) = struct();
        results_baseline.(sample_name).data = data;
        results_baseline.(sample_name).res = res_baseline;
        results_baseline.(sample_name).irfs = irfs_baseline;
        results_baseline.(sample_name).shocknames = {'Pure MP', 'CBI'};
    end
    
    %% --------------------------------------------------------------------
    % POOR MAN'S: Cholesky on (PMNEG, PMPOS) proxies
    % --------------------------------------------------------------------
    if RUN_POORMAN
        fprintf('\n--- POOR MAN''S: Cholesky on proxies ---\n');
        
        % Prepare data structure
        data_pm = struct();
        
        % Use Poor Man's proxies as i.i.d. instruments
        data_pm.y = [sample_data.pmneg_ed2_sp500_tc, ...  % m1: Pure MP proxy
                     sample_data.pmpos_ed2_sp500_tc, ...  % m2: CBI proxy
                     sample_data.gs1, ...                  % y1: 1-year rate
                     sample_data.logcpi, ...               % y2: log CPI
                     sample_data.logip, ...                % y3: log IP
                     sample_data.ebp];                     % y4: EBP
        
        data_pm.names = {'pmneg_ed2_sp500_tc', 'pmpos_ed2_sp500_tc', 'gs1', 'logcpi', 'logip', 'ebp'};
        data_pm.Nm = 2;
        data_pm.w = ones(size(data_pm.y, 1), 1);
        data_pm.time = ym2t(sample_data.year, sample_data.month);
        
        % Prior settings
        prior.Nm = data_pm.Nm;
        prior.minnesota.mvector = [0, 0, 1, 0, 0, 1];
        
        % Diagnostics
        fprintf('  PMNEG NaNs: %d, Zeros: %d (%.1f%% valid)\n', ...
            sum(isnan(data_pm.y(:,1))), sum(data_pm.y(:,1)==0), ...
            100*mean(~isnan(data_pm.y(:,1)) & data_pm.y(:,1)~=0));
        fprintf('  PMPOS NaNs: %d, Zeros: %d (%.1f%% valid)\n', ...
            sum(isnan(data_pm.y(:,2))), sum(data_pm.y(:,2)==0), ...
            100*mean(~isnan(data_pm.y(:,2)) & data_pm.y(:,2)~=0));
        
        % Orthogonality check (JK 2020 emphasize distinct structural shocks)
        valid_idx = ~isnan(data_pm.y(:,1)) & ~isnan(data_pm.y(:,2)) & ...
                    data_pm.y(:,1)~=0 & data_pm.y(:,2)~=0;
        if sum(valid_idx) > 1
            corr_mp_cbi = corr(data_pm.y(valid_idx,1), data_pm.y(valid_idx,2));
            fprintf('  Correlation(PMNEG, PMPOS): %.3f (should be small for distinct shocks)\n', corr_mp_cbi);
        end
        
        % Estimate VAR
        fprintf('  Estimating VAR... ');
        res_pm = VAR_withiid1kf(data_pm, prior, gssettings);
        fprintf('Done!\n');
        
        % Cholesky identification (JK's "Poor Man's" approach)
        fprintf('  Computing Cholesky IRFs... ');
        N_pm = length(data_pm.names);
        irfs_pm = NaN(N_pm, N_pm, MAlags, gssettings.ndraws);
        
        for i = 1:gssettings.ndraws
            betadraw = res_pm.beta_draws(1:end-size(data_pm.w,2),:,i);
            sigmadraw = res_pm.sigma_draws(:,:,i);
            response = impulsdtrf(reshape(betadraw', N_pm, N_pm, prior.lags), chol(sigmadraw), MAlags);
            irfs_pm(:,:,:,i) = response;
        end
        fprintf('Done!\n');
        
        % Store results
        results_pm.(sample_name) = struct();
        results_pm.(sample_name).data = data_pm;
        results_pm.(sample_name).res = res_pm;
        results_pm.(sample_name).irfs = irfs_pm;
        results_pm.(sample_name).shocknames = {'Pure MP', 'CBI'};
    end
end

%% ========================================================================
% SAVE RESULTS
% =========================================================================

fprintf('\n========================================\n');
fprintf('Saving results...\n');
fprintf('========================================\n');

if RUN_BASELINE
    save(fullfile(output_dir, 'JK_Baseline_Results.mat'), 'results_baseline', '-v7.3');
    fprintf('  Baseline results saved to: %s\n', fullfile(output_dir, 'JK_Baseline_Results.mat'));
end

if RUN_POORMAN
    save(fullfile(output_dir, 'JK_PoorMan_Results.mat'), 'results_pm', '-v7.3');
    fprintf('  Poor Man''s results saved to: %s\n', fullfile(output_dir, 'JK_PoorMan_Results.mat'));
end

fprintf('\n✓ All estimations complete!\n');
fprintf('\nNext steps:\n');
fprintf('  1. Run plot_JK_results.m to generate figures\n');
fprintf('  2. Compare Baseline vs Poor Man''s decompositions\n');
fprintf('  3. Compare across sample periods\n\n');
