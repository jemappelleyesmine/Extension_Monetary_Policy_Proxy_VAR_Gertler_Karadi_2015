% =========================================================================
% PLOT JK (2020) DECOMPOSITION RESULTS (IMPROVED VERSION)
% Generates comparison figures for Pure MP vs CBI shocks
% 
% IMPROVEMENTS FROM REVIEWER FEEDBACK:
%   1. Explicit variable indexing (no hard-coded positions)
%   2. Clear normalization statement
%   3. Diagnostic output (proxy statistics)
%   4. Consistent confidence band terminology
%   5. Improved figure layouts (JK-style comparison)
% =========================================================================

clear all; close all;

%% ========================================================================
% CONFIGURATION
% =========================================================================

% Output directory
output_dir = 'output';
if ~exist(output_dir, 'dir')
    error('Output directory not found. Run run_JK_decomposition.m first.');
end

% Which results to plot
PLOT_BASELINE = true;   % Baseline sign restrictions
PLOT_POORMAN = true;    % Poor Man's (Cholesky)
PLOT_COMPARISON = true; % Side-by-side comparison

% Plot settings
qtoplot = [0.5 0.16 0.84];  % Median and 68% credible sets
horiz_plot = 48;  % Months to plot

% NORMALIZATION NOTE:
% All shocks are normalized to 1 standard deviation (unit variance).
% This follows JK (2020) and allows direct comparison across methodologies.
fprintf('\n========================================\n');
fprintf('NORMALIZATION: All responses are to 1 SD shocks\n');
fprintf('========================================\n\n');

%% ========================================================================
% LOAD RESULTS AND CREATE VARIABLE INDEX MAP
% =========================================================================

if PLOT_BASELINE
    fprintf('Loading baseline results...\n');
    load(fullfile(output_dir, 'JK_Baseline_Results.mat'), 'results_baseline');
    samples_baseline = fieldnames(results_baseline);
    
    % Get variable names from first sample
    sample1 = samples_baseline{1};
    varnames = results_baseline.(sample1).data.names;
    fprintf('  Variables in order: %s\n', strjoin(varnames, ', '));
end

if PLOT_POORMAN
    fprintf('Loading Poor Man''s results...\n');
    load(fullfile(output_dir, 'JK_PoorMan_Results.mat'), 'results_pm');
    samples_pm = fieldnames(results_pm);
    
    if ~exist('varnames', 'var')
        sample1 = samples_pm{1};
        varnames = results_pm.(sample1).data.names;
        fprintf('  Variables in order: %s\n', strjoin(varnames, ', '));
    end
end

%% ========================================================================
% CREATE EXPLICIT VARIABLE INDEX MAP
% =========================================================================
% This ensures we always plot the correct variables regardless of ordering

fprintf('\nCreating variable index map...\n');

% Define the variables to plot (macro variables only, not instruments)
vars_to_plot = {'gs1', 'logcpi', 'logip', 'ebp'};
var_labels = {'1-Year Rate', 'CPI', 'Industrial Production', 'EBP'};

% Create index map
idx = struct();
for i = 1:length(vars_to_plot)
    varname = vars_to_plot{i};
    pos = find(strcmpi(varnames, varname));
    
    if isempty(pos)
        error('Variable %s not found in results!', varname);
    end
    
    idx.(varname) = pos;
    fprintf('  %s -> position %d\n', varname, pos);
end

% Shock names
shock_names = {'Pure MP', 'CBI'};

%% ========================================================================
% PRINT PROXY DIAGNOSTICS
% =========================================================================

if PLOT_POORMAN
    fprintf('\n========================================\n');
    fprintf('PROXY DIAGNOSTICS (Poor Man''s Method)\n');
    fprintf('========================================\n');
    
    for s = 1:length(samples_pm)
        sample_name = samples_pm{s};
        data_pm = results_pm.(sample_name).data;
        
        % Extract proxies (first two variables)
        pmneg = data_pm.y(:,1);
        pmpos = data_pm.y(:,2);
        
        % Compute statistics
        pmneg_nonzero = mean(~isnan(pmneg) & pmneg ~= 0) * 100;
        pmpos_nonzero = mean(~isnan(pmpos) & pmpos ~= 0) * 100;
        
        % Correlation (only valid observations)
        valid_idx = ~isnan(pmneg) & ~isnan(pmpos) & pmneg ~= 0 & pmpos ~= 0;
        if sum(valid_idx) > 1
            proxy_corr = corr(pmneg(valid_idx), pmpos(valid_idx));
        else
            proxy_corr = NaN;
        end
        
        fprintf('\n%s:\n', sample_name);
        fprintf('  MP proxy non-zero: %.1f%%\n', pmneg_nonzero);
        fprintf('  CBI proxy non-zero: %.1f%%\n', pmpos_nonzero);
        fprintf('  Correlation(MP, CBI): %.3f\n', proxy_corr);
    end
    fprintf('\n');
end

%% ========================================================================
% HELPER FUNCTION: Extract IRF quantiles
% =========================================================================

function [irf_median, irf_low, irf_high] = extract_irfs(irfs_draws, var_idx, shock_idx, quantiles)
    % Extract IRFs for specific variable and shock
    % irfs_draws: N x N x H x ndraws
    % Returns: 1 x H vectors
    
    irf_data = squeeze(irfs_draws(var_idx, shock_idx, :, :));  % H x ndraws
    
    irf_median = quantile(irf_data, quantiles(1), 2)';
    irf_low = quantile(irf_data, quantiles(2), 2)';
    irf_high = quantile(irf_data, quantiles(3), 2)';
end

%% ========================================================================
% PLOT BASELINE RESULTS (One figure per sample)
% =========================================================================

if PLOT_BASELINE
    fprintf('\n========================================\n');
    fprintf('GENERATING BASELINE FIGURES\n');
    fprintf('========================================\n');
    
    for s = 1:length(samples_baseline)
        sample_name = samples_baseline{s};
        
        fprintf('\nPlotting baseline: %s\n', sample_name);
        
        irfs = results_baseline.(sample_name).irfs;
        
        % Create figure
        figure('Position', [100 100 1200 800]);
        sgtitle(sprintf('JK Baseline: %s', sample_name), 'FontSize', 14, 'FontWeight', 'bold');
        
        % Loop over variables using explicit index
        for v = 1:length(vars_to_plot)
            varname = vars_to_plot{v};
            var_idx = idx.(varname);
            
            for sh = 1:2  % MP and CBI shocks
                subplot(length(vars_to_plot), 2, (v-1)*2 + sh);
                
                % Extract IRFs
                [irf_med, irf_low, irf_high] = extract_irfs(irfs, var_idx, sh, qtoplot);
                
                % Plot
                x = 0:horiz_plot-1;
                fill([x fliplr(x)], [irf_high(1:horiz_plot) fliplr(irf_low(1:horiz_plot))], ...
                    [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                hold on;
                plot(x, irf_med(1:horiz_plot), 'b-', 'LineWidth', 2);
                plot(x, zeros(1,horiz_plot), 'k--', 'LineWidth', 0.5);
                hold off;
                
                % Format
                xlabel('Months');
                if sh == 1
                    ylabel(var_labels{v});
                end
                title(shock_names{sh});
                grid on;
                xlim([0 horiz_plot-1]);
            end
        end
        
        % Add note about normalization and bands
        annotation('textbox', [0.01 0.01 0.98 0.03], ...
            'String', 'Note: Responses to 1 SD shock. Shaded areas show 68% posterior credible sets.', ...
            'EdgeColor', 'none', 'FontSize', 8, 'HorizontalAlignment', 'center');
        
        % Save figure
        filename = fullfile(output_dir, sprintf('JK_Baseline_%s.png', sample_name));
        saveas(gcf, filename);
        fprintf('  Saved: %s\n', filename);
    end
end

%% ========================================================================
% PLOT POOR MAN'S RESULTS (One figure per sample)
% =========================================================================

if PLOT_POORMAN
    fprintf('\n========================================\n');
    fprintf('GENERATING POOR MAN''S FIGURES\n');
    fprintf('========================================\n');
    
    for s = 1:length(samples_pm)
        sample_name = samples_pm{s};
        
        fprintf('\nPlotting Poor Man''s: %s\n', sample_name);
        
        irfs = results_pm.(sample_name).irfs;
        
        % Create figure
        figure('Position', [100 100 1200 800]);
        sgtitle(sprintf('JK Poor Man: %s', sample_name), 'FontSize', 14, 'FontWeight', 'bold');
        
        % Loop over variables using explicit index
        for v = 1:length(vars_to_plot)
            varname = vars_to_plot{v};
            var_idx = idx.(varname);
            
            for sh = 1:2  % MP and CBI shocks
                subplot(length(vars_to_plot), 2, (v-1)*2 + sh);
                
                % Extract IRFs
                [irf_med, irf_low, irf_high] = extract_irfs(irfs, var_idx, sh, qtoplot);
                
                % Plot
                x = 0:horiz_plot-1;
                fill([x fliplr(x)], [irf_high(1:horiz_plot) fliplr(irf_low(1:horiz_plot))], ...
                    [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                hold on;
                plot(x, irf_med(1:horiz_plot), 'r-', 'LineWidth', 2);
                plot(x, zeros(1,horiz_plot), 'k--', 'LineWidth', 0.5);
                hold off;
                
                % Format
                xlabel('Months');
                if sh == 1
                    ylabel(var_labels{v});
                end
                title(shock_names{sh});
                grid on;
                xlim([0 horiz_plot-1]);
            end
        end
        
        % Add note
        annotation('textbox', [0.01 0.01 0.98 0.03], ...
            'String', 'Note: Responses to 1 SD shock. Shaded areas show 68% bootstrap percentile bands.', ...
            'EdgeColor', 'none', 'FontSize', 8, 'HorizontalAlignment', 'center');
        
        % Save figure
        filename = fullfile(output_dir, sprintf('JK_PoorMan_%s.png', sample_name));
        saveas(gcf, filename);
        fprintf('  Saved: %s\n', filename);
    end
end

%% ========================================================================
% COMPARISON PLOTS: Baseline vs Poor Man's (for each sample)
% =========================================================================

if PLOT_BASELINE && PLOT_POORMAN && PLOT_COMPARISON
    fprintf('\n========================================\n');
    fprintf('GENERATING COMPARISON FIGURES\n');
    fprintf('========================================\n');
    
    % Compare all samples that exist in both results
    samples_common = intersect(samples_baseline, samples_pm);
    
    for s = 1:length(samples_common)
        sample_name = samples_common{s};
        
        fprintf('\nPlotting comparison: %s\n', sample_name);
        
        figure('Position', [100 100 1400 900]);
        sgtitle(sprintf('Comparison: Baseline vs Poor Man (%s)', sample_name), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        % Loop over variables
        for v = 1:length(vars_to_plot)
            varname = vars_to_plot{v};
            var_idx = idx.(varname);
            
            for sh = 1:2
                subplot(length(vars_to_plot), 2, (v-1)*2 + sh);
                
                % Baseline IRFs
                [irf_base_med, irf_base_low, irf_base_high] = ...
                    extract_irfs(results_baseline.(sample_name).irfs, var_idx, sh, qtoplot);
                
                % Poor Man's IRFs
                [irf_pm_med, irf_pm_low, irf_pm_high] = ...
                    extract_irfs(results_pm.(sample_name).irfs, var_idx, sh, qtoplot);
                
                % Plot
                x = 0:horiz_plot-1;
                
                % Baseline (blue)
                fill([x fliplr(x)], [irf_base_high(1:horiz_plot) fliplr(irf_base_low(1:horiz_plot))], ...
                    [0.7 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                hold on;
                plot(x, irf_base_med(1:horiz_plot), 'b-', 'LineWidth', 2);
                
                % Poor Man's (red)
                fill([x fliplr(x)], [irf_pm_high(1:horiz_plot) fliplr(irf_pm_low(1:horiz_plot))], ...
                    [1 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                plot(x, irf_pm_med(1:horiz_plot), 'r-', 'LineWidth', 2);
                
                plot(x, zeros(1,horiz_plot), 'k--', 'LineWidth', 0.5);
                hold off;
                
                % Format
                xlabel('Months');
                if sh == 1
                    ylabel(var_labels{v});
                end
                title(shock_names{sh});
                if v == 1 && sh == 1
                    legend({'Baseline 68%', 'Baseline', 'Poor Man 68%', 'Poor Man'}, ...
                        'Location', 'best', 'FontSize', 8);
                end
                grid on;
                xlim([0 horiz_plot-1]);
            end
        end
        
        % Add note
        annotation('textbox', [0.01 0.01 0.98 0.03], ...
            'String', 'Note: All responses to 1 SD shocks. Blue = Baseline (posterior), Red = Poor Man (bootstrap).', ...
            'EdgeColor', 'none', 'FontSize', 8, 'HorizontalAlignment', 'center');
        
        % Save
        filename = fullfile(output_dir, sprintf('JK_Comparison_%s.png', sample_name));
        saveas(gcf, filename);
        fprintf('  Saved: %s\n', filename);
    end
end

%% ========================================================================
% ALTERNATIVE LAYOUT: JK-STYLE THREE-PANEL COMPARISON
% (All three methodologies on same axes)
% =========================================================================

if PLOT_BASELINE && PLOT_POORMAN
    fprintf('\n========================================\n');
    fprintf('GENERATING JK-STYLE OVERLAY FIGURES\n');
    fprintf('========================================\n');
    
    samples_common = intersect(samples_baseline, samples_pm);
    
    for s = 1:length(samples_common)
        sample_name = samples_common{s};
        
        fprintf('\nPlotting JK-style overlay: %s\n', sample_name);
        
        figure('Position', [100 100 1400 900]);
        sgtitle(sprintf('MP vs CBI Shocks (%s)', sample_name), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        % Loop over variables
        for v = 1:length(vars_to_plot)
            varname = vars_to_plot{v};
            var_idx = idx.(varname);
            
            for sh = 1:2
                subplot(length(vars_to_plot), 2, (v-1)*2 + sh);
                
                % Extract both Baseline and Poor Man's
                [irf_base, ~, ~] = extract_irfs(results_baseline.(sample_name).irfs, var_idx, sh, qtoplot);
                [irf_pm, ~, ~] = extract_irfs(results_pm.(sample_name).irfs, var_idx, sh, qtoplot);
                
                x = 0:horiz_plot-1;
                
                % Plot both lines
                plot(x, irf_base(1:horiz_plot), 'b-', 'LineWidth', 2.5);
                hold on;
                plot(x, irf_pm(1:horiz_plot), 'r--', 'LineWidth', 2);
                plot(x, zeros(1,horiz_plot), 'k--', 'LineWidth', 0.5);
                hold off;
                
                xlabel('Months');
                if sh == 1
                    ylabel(var_labels{v});
                end
                title(shock_names{sh});
                if v == 1 && sh == 1
                    legend({'Baseline (Sign)', 'Poor Man (Chol)'}, 'Location', 'best');
                end
                grid on;
                xlim([0 horiz_plot-1]);
            end
        end
        
        % Add note
        annotation('textbox', [0.01 0.01 0.98 0.03], ...
            'String', 'Note: Median IRFs to 1 SD shocks. Blue = Baseline sign restrictions, Red = Poor Man Cholesky.', ...
            'EdgeColor', 'none', 'FontSize', 8, 'HorizontalAlignment', 'center');
        
        filename = fullfile(output_dir, sprintf('JK_Overlay_%s.png', sample_name));
        saveas(gcf, filename);
        fprintf('  Saved: %s\n', filename);
    end
end

fprintf('\n========================================\n');
fprintf('✓ All figures generated!\n');
fprintf('========================================\n\n');

fprintf('FIGURE SUMMARY:\n');
fprintf('  - Individual baseline figures (per sample)\n');
fprintf('  - Individual Poor Man figures (per sample)\n');
fprintf('  - Side-by-side comparisons (per sample)\n');
fprintf('  - JK-style overlay figures (per sample)\n\n');

fprintf('NEXT STEPS:\n');
fprintf('  1. Review figures in output/ folder\n');
fprintf('  2. Choose which figures to include in report\n');
fprintf('  3. Use METHODOLOGY_NOTES.md for caption templates\n\n');
