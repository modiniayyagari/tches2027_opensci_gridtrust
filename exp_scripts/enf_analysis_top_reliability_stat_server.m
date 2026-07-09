%% ENF Server-Room Robustness Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Evaluates whether the ENF signature captured on a DC-powered server inside a
%   server room remains correlated with the AC mains reference. The analysis spans
%   two sensor placement configurations (SBOX and SPSU) and two days-of-week
%   (MON and WED) to test robustness across sensor placement and temporal conditions.
%
% Analysis Steps:
%   1. Scan exp_inputs/SRV_L/ for matching mains_pow_trace_ac.wav and
%      fpga_em_trace_dc.wav pairs under the SBOX/ and SPSU/ sub-trees.
%   2. Compute the Pearson correlation (r) for each trace pair using the
%      shared ENF extraction pipeline (proc_enf_analysis).
%   3. Apply the Fisher z = atanh(r) transform for variance stabilization.
%   4. Report descriptive statistics (N, mean, SD, 95% CI) for r and z.
%   5. Run one-way repeated-measures ANOVA for each factor:
%        - Day effect:  MON vs WED, with Time as the repeated-measures
%                       subject (correlations averaged over Site before pivoting).
%        - Site effect: SBOX vs SPSU, with Time as the repeated-measures
%                       subject (correlations averaged over Day before pivoting).
%   6. Export descriptive box plots, histograms, z-distribution curves, and
%      dot-and-whisker mean-CI plots to exp_results/srv_room_robustness/.
%
% Inputs:
%   exp_inputs/SRV_L/ - Artifact dataset with SBOX/ and SPSU/ sub-trees,
%                       each holding MON/ and WED/ day folders with T01-T05
%                       trial leaf folders containing paired .wav trace files.
%
% Outputs:
%   exp_results/srv_room_robustness/  - All plot files (.png, .svg, .fig)
%   exp_logs/<script_name>_log.txt    - Full text log of all printed output
%
% Dependencies:
%   proc_enf_analysis - Pre-compiled ENF extraction and correlation function.
%                       Must be on the MATLAB path before running this script.
%
% Usage:
%   Run from the exp_scripts/ directory. The artifact data root is resolved
%   automatically relative to this script file's location.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% 1) Script Configuration
% Configure all analysis parameters in the struct below before running.
% Parameters cover dataset paths, ENF extraction settings, and output flags.
config = struct();

% Resolve the artifact data root relative to this script location.
script_dir    = fileparts(mfilename('fullpath'));
artifact_root = fileparts(script_dir);
config.artifact_root = artifact_root;
config.baseDir = fullfile(artifact_root, 'exp_inputs', 'SRV_L');

% WAV base names used in each SRV_L trial folder.
config.file_1_name_base = 'mains_pow_trace_ac'; % reference
config.file_2_name_base = 'fpga_em_trace_dc';   % sensed

% Expected dataset labels.
config.allowed_sites = {'SBOX','SPSU'};
config.allowed_days  = {'MON','WED'};
config.allowed_times = {'T01','T02','T03','T04','T05'};

% Reporting flags.
config.SCR_CFG_DISP_RAW_PAIR_WISE_CORR_VALS = false;
config.SCR_CFG_DISP_DESCRIPTIVE_STATS_RES   = true;
config.SCR_CFG_LIST_SUBJECT_KEYS            = false;

% ENF extraction parameters
config.nominal_freq_arr = [50 60];  % use 60 Hz
config.frame_size       = 8000;     % ms
config.nfft             = 2^15;
config.overlap_size     = 0;        % ms

nominal_freq = config.nominal_freq_arr(2); % 60 Hz
harmonics    = (1:7) * nominal_freq;

config.harmonics_arr_1  = harmonics;
config.trace_1_freq_est_method = 1;
config.trace_1_est_freq = nominal_freq;
config.trace_1_freq_est_spec_comb_harmonics = harmonics;
config.trace_1_plot_title = config.file_1_name_base;

config.harmonics_arr_2  = harmonics;
config.trace_2_freq_est_method = 1;
config.trace_2_est_freq = nominal_freq;
config.trace_2_freq_est_spec_comb_harmonics = harmonics;
config.trace_2_plot_title = config.file_2_name_base;

%% 2) Setup: Logging and Output Directories
% Create the output and log directories if they do not exist, open the log
% file for writing, and register a cleanup handler to close it on exit.
[fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders_srv();
config.results_dir = results_dir;

fprintf('Starting ENF Server-Room Robustness Analysis. Log: %s\n', log_filename);
fprintf(fileID, 'Starting ENF Server-Room Robustness Analysis. Log: %s\n\n', log_filename);

%% 3) ENF Correlation Batch: Scan and Compute
% Recursively scan exp_inputs/SRV_L/ for trace pairs, call proc_enf_analysis
% on each pair, and collect results into a table (Correlation, Site, Day, Time).
tbl = run_serverroom_scan(config, fileID);  % Correlation, FilePath, Site, Day, Time

if isempty(tbl)
    fprintf('No data found. Exiting.\n');
    fprintf(fileID, 'No data found. Exiting.\n');
    return;
end

all_r = tbl.Correlation;

%% 4) Descriptive Statistics and Summary Plots
% Compute r-scale and Fisher z-scale descriptive statistics for all collected
% correlations. Generate box plots, histograms, and z-score distribution
% figures and save them to the results directory.
[r_stats, z_data, z_stats] = calculate_basic_stats(all_r);

if config.SCR_CFG_DISP_DESCRIPTIVE_STATS_RES
    print_descriptive_stats_results(r_stats, z_stats, fileID);
    try
        plot_box_simple(all_r, r_stats, results_dir, fileID, 'r');
        plot_hist_simple(all_r, r_stats, results_dir, fileID, 'r');
        plot_box_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_hist_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_z_score_distributions(z_data, z_stats, results_dir, fileID);
    catch MEp
        fprintf('Plot error: %s\n', MEp.message);
        fprintf(fileID, 'Plot error: %s\n', MEp.message);
    end
else
    fprintf('\nSTEP 1 (Descriptive) skipped by config.\n');
    fprintf(fileID, '\nSTEP 1 (Descriptive) skipped by config.\n');
end

if config.SCR_CFG_DISP_RAW_PAIR_WISE_CORR_VALS
    print_raw_correlation_data_noID(tbl, all_r, fileID);
end

%% 5) Inferential Test: Day Factor (MON vs WED)
% Run a one-way repeated-measures ANOVA comparing MON and WED correlations.
% Correlations are averaged over Site within each Time, so Time is the
% repeated-measures subject. Reports F, df, p, p(GG), and eta_p^2.
fprintf('\n[DAY] 1-way RM-ANOVA (MON vs WED), subjects = Time\n');
fprintf(fileID, '\n[DAY] 1-way RM-ANOVA (MON vs WED), subjects = Time\n');
run_oneway_rm_day(tbl, config, fileID);

%% 6) Inferential Test: Site Factor (SBOX vs SPSU)
% Run a one-way repeated-measures ANOVA comparing SBOX and SPSU correlations.
% Correlations are averaged over Day within each Time, so Time is the
% repeated-measures subject. Reports F, df, p, p(GG), and eta_p^2.
fprintf('\n[SITE] 1-way RM-ANOVA (SBOX vs SPSU), subjects = Time\n');
fprintf(fileID, '\n[SITE] 1-way RM-ANOVA (SBOX vs SPSU), subjects = Time\n');
run_oneway_rm_site(tbl, config, fileID);

fprintf('\nDone. Results in: %s\n', results_dir);
fprintf(fileID, '\nDone. Results in: %s\n', results_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders_srv()
    log_dir     = 'exp_logs';
    results_dir = 'exp_results/srv_room_robustness';
    if ~exist(log_dir,'dir'),     mkdir(log_dir);     end
    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    [~, script_name] = fileparts(mfilename('fullpath'));
    if isempty(script_name), script_name = 'enf_srv_room_robustness'; end
    log_filename = fullfile(log_dir, [script_name '_log.txt']);
    if exist(log_filename,'file'), delete(log_filename); end
    fileID = fopen(log_filename,'w');
    if fileID == -1, error('Could not open %s', log_filename); end
    cleanupObj = onCleanup(@() fclose(fileID));
end

% run_serverroom_scan - Scan SRV_L folders and compute pairwise ENF correlations.
%
%   Recursively searches config.baseDir for all fpga_em_trace_dc.wav files.
%   For each sensed file, the co-located mains_pow_trace_ac.wav reference is
%   matched and proc_enf_analysis is called to compute the Pearson correlation.
%   Site, Day, and Time labels are parsed from path tokens using regex.
%
%   Inputs:
%     config  - Configuration struct (baseDir, file names, ENF parameters, flags)
%     fileID  - Log file identifier for tee-printing to the log file
%
%   Outputs:
%     tbl - Table with columns: Correlation, FilePath, Site, Day, Time
%
function tbl = run_serverroom_scan(config, fileID)
    varTypes = {'double','string','string','string','string'};
    varNames = {'Correlation','FilePath','Site','Day','Time'};

    scan_label = artifact_relative_path(config.baseDir, config.artifact_root);
    fprintf('INFO: Scanning %s\n', scan_label);
    fprintf(fileID, 'INFO: Scanning %s\n', scan_label);

    sensed_files = dir(fullfile(config.baseDir, '**', [config.file_2_name_base, '.wav']));
    fprintf('INFO: Found %d sensed files.\n', numel(sensed_files));
    fprintf(fileID, 'INFO: Found %d sensed files.\n', numel(sensed_files));

    tbl = table('Size',[numel(sensed_files),numel(varNames)], ...
                'VariableTypes',varTypes, ...
                'VariableNames',varNames);
    rowCount = 0;

    % Case-insensitive token finder
    hasTok = @(p,tok) ~isempty(regexpi(p, ['(^|[\\/_.-])' tok '([\\/_.-]|$)'], 'once'));

    for k = 1:numel(sensed_files)
        thisFolder = sensed_files(k).folder;

        % Parse Site / Day / Time from path
        site = '';
        for i = 1:numel(config.allowed_sites)
            if hasTok(thisFolder, config.allowed_sites{i})
                site = upper(config.allowed_sites{i});
                break;
            end
        end

        day = '';
        for i = 1:numel(config.allowed_days)
            if hasTok(thisFolder, config.allowed_days{i})
                day = upper(config.allowed_days{i});
                break;
            end
        end

        time = '';
        for i = 1:numel(config.allowed_times)
            if hasTok(thisFolder, config.allowed_times{i})
                time = upper(config.allowed_times{i});
                break;
            end
        end

        if any(cellfun(@isempty, {site, day, time}))
            continue;
        end

        f2 = fullfile(thisFolder, [config.file_2_name_base, '.wav']); % sensed
        f1 = fullfile(thisFolder, [config.file_1_name_base, '.wav']); % reference

        if ~isfile(f1)
            f2_label = artifact_relative_path(f2, config.artifact_root);
            fprintf('WARN: Missing reference WAV for %s\n', f2_label);
            fprintf(fileID, 'WARN: Missing reference WAV for %s\n', f2_label);
            continue;
        end

        try
            % Compute ENF-based correlation using external pipeline
            r = proc_enf_analysis( ...
                    f1, f2, ...
                    config.nfft, config.frame_size, config.overlap_size, ...
                    config.harmonics_arr_1, config.nominal_freq_arr(2), ...
                    config.harmonics_arr_2, config.nominal_freq_arr(2), ...
                    config.trace_1_freq_est_method, config.trace_1_est_freq, config.trace_1_freq_est_spec_comb_harmonics, ...
                    config.trace_2_freq_est_method, config.trace_2_est_freq, config.trace_2_freq_est_spec_comb_harmonics, ...
                    config.trace_1_plot_title, config.trace_2_plot_title, false);

            rowCount = rowCount + 1;
            tbl(rowCount,:) = {r, artifact_relative_path(f2, config.artifact_root), ...
                               string(site), string(day), string(time)};

            fprintf('INFO: %3d/%d  %s | %s | %s   r=%.6f\n', ...
                    k, numel(sensed_files), site, day, time, r);
            fprintf(fileID, 'INFO: %3d/%d  %s | %s | %s   r=%.6f\n', ...
                    k, numel(sensed_files), site, day, time, r);
        catch ME
            fprintf('ERROR: %s\n', ME.message);
            fprintf(fileID, 'ERROR: %s\n', ME.message);
        end
    end

    tbl = tbl(1:rowCount, :);
    fprintf('\nINFO: Batch complete. N=%d\n', height(tbl));
    fprintf(fileID, '\nINFO: Batch complete. N=%d\n', height(tbl));
end

function rel_path = artifact_relative_path(path_in, artifact_root)
    path_str = strrep(char(path_in), '\', '/');
    root_str = strrep(char(artifact_root), '\', '/');
    if numel(path_str) > numel(root_str) && ...
            strcmpi(path_str(1:numel(root_str)), root_str) && ...
            path_str(numel(root_str) + 1) == '/'
        rel_path = string(path_str(numel(root_str) + 2:end));
    elseif strcmpi(path_str, root_str)
        rel_path = ".";
    else
        rel_path = string(path_str);
    end
end

% calculate_basic_stats - Compute r-scale and Fisher z-scale descriptive statistics.
%
%   Clips r values to the open interval (-1, 1) before applying atanh() so
%   that boundary values do not produce Inf in the z-transform.
%
%   Inputs:
%     all_r - Numeric vector of Pearson correlation coefficients
%
%   Outputs:
%     r_stats - Struct of descriptive statistics for the raw r values
%     z_data  - Fisher z-transformed values (z = atanh(r), clipped)
%     z_stats - Struct of descriptive statistics for the z values
%
function [r_stats, z_data, z_stats] = calculate_basic_stats(all_r)
    r_stats = sub_calc_stats(all_r);
    r_safe = all_r;
    r_safe(r_safe >=  1) = 0.9999999;
    r_safe(r_safe <= -1) = -0.9999999;
    z_data  = atanh(r_safe);
    z_stats = sub_calc_stats(z_data);
end

% sub_calc_stats - Compute a full set of descriptive statistics for a data vector.
%
%   Returns count, central tendency (mean, median), spread (std, variance,
%   skewness, kurtosis, SEM), and the 95% CI margin using the t-critical value
%   for the sample size.
%
%   Inputs:
%     x - Numeric vector of data values (may be empty)
%
%   Outputs:
%     stats - Struct with fields: N, mean, median, std, min, max, var,
%             skew, kurt, sem, ci95_margin
%
function stats = sub_calc_stats(x)
    if ~isempty(x)
        stats.N      = numel(x);
        stats.mean   = mean(x);
        stats.median = median(x);
        stats.std    = std(x);
        stats.min    = min(x);
        stats.max    = max(x);
        stats.var    = var(x);
        stats.skew   = skewness(x);
        stats.kurt   = kurtosis(x);
        stats.sem    = stats.std / sqrt(stats.N);
        df = max(stats.N - 1, 1);
        tcrit = tinv(0.975, df);
        stats.ci95_margin = tcrit * stats.sem;
    else
        stats = struct('N',0,'mean',NaN,'median',NaN,'std',NaN,'min',NaN,'max',NaN, ...
                       'var',NaN,'skew',NaN,'kurt',NaN,'sem',NaN,'ci95_margin',NaN);
    end
end

% print_descriptive_stats_results - Print r-scale and z-scale descriptive stats.
%
%   Writes STEP 1A (raw r) and STEP 1B (Fisher z) descriptive tables to both
%   the MATLAB console and the log file.
%
%   Inputs:
%     r_stats - Stats struct for raw Pearson r values (from sub_calc_stats)
%     z_stats - Stats struct for Fisher z values (from sub_calc_stats)
%     fileID  - Log file identifier for tee-printing
%
function print_descriptive_stats_results(r_stats, z_stats, fileID)
    fprintf('\nSTEP 1A: Descriptive Results (r)\n');
    fprintf(fileID, '\nSTEP 1A: Descriptive Results (r)\n');
    print_stats_table(r_stats, fileID, 'All Correlations (r)');

    fprintf('\nSTEP 1B: Descriptive Results (z)\n');
    fprintf(fileID, '\nSTEP 1B: Descriptive Results (z)\n');
    print_stats_table(z_stats, fileID, 'All Correlations (z'')');

    fprintf('STEP 1 complete.\n');
    fprintf(fileID, 'STEP 1 complete.\n');
end

% print_stats_table - Print a single-column descriptive statistics table.
%
%   Writes N, mean, median, SD, variance, skewness, kurtosis, SEM, 95% CI,
%   min, and max to both the console and the log file. When the label string
%   contains "z'", the 95% CI is back-transformed to r-scale via tanh and
%   printed as a second CI row.
%
%   Inputs:
%     stats  - Stats struct (from sub_calc_stats)
%     fileID - Log file identifier for tee-printing
%     label  - Column header label string (e.g., 'All Correlations (r)')
%
function print_stats_table(stats, fileID, label)
    fprintf('%-12s | %20s\n', 'Statistic', label);
    fprintf(fileID, '%-12s | %20s\n', 'Statistic', label);

    fprintf('%-12s | %20d\n', 'Count (N)', stats.N);
    fprintf(fileID, '%-12s | %20d\n', 'Count (N)', stats.N);

    fprintf('%-12s | %20.4f\n', 'Mean', stats.mean);
    fprintf(fileID, '%-12s | %20.4f\n', 'Mean', stats.mean);

    fprintf('%-12s | %20.4f\n', 'Median', stats.median);
    fprintf(fileID, '%-12s | %20.4f\n', 'Median', stats.median);

    fprintf('%-12s | %20.4f\n', 'Std. Dev.', stats.std);
    fprintf(fileID, '%-12s | %20.4f\n', 'Std. Dev.', stats.std);

    fprintf('%-12s | %20.4f\n', 'Variance', stats.var);
    fprintf(fileID, '%-12s | %20.4f\n', 'Variance', stats.var);

    fprintf('%-12s | %20.4f\n', 'Skewness', stats.skew);
    fprintf(fileID, '%-12s | %20.4f\n', 'Skewness', stats.skew);

    fprintf('%-12s | %20.4f\n', 'Kurtosis', stats.kurt);
    fprintf(fileID, '%-12s | %20.4f\n', 'Kurtosis', stats.kurt);

    fprintf('%-12s | %20.4f\n', 'Std. Error', stats.sem);
    fprintf(fileID, '%-12s | %20.4f\n', 'Std. Error', stats.sem);

    ci_str = sprintf('[%.4f, %.4f]', ...
        stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin);
    fprintf('%-12s | %20s\n', '95% CI', ci_str);
    fprintf(fileID, '%-12s | %20s\n', '95% CI', ci_str);

    if contains(label,"z'")
        ci_r = tanh([stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin]);
        ci_r_str = sprintf('[%.4f, %.4f]', ci_r(1), ci_r(2));
        fprintf('%-12s | %20s\n', '  (95% CI r)', ci_r_str);
        fprintf(fileID, '%-12s | %20s\n', '  (95% CI r)', ci_r_str);
    end

    fprintf('%-12s | %20.4f\n', 'Min', stats.min);
    fprintf(fileID, '%-12s | %20.4f\n', 'Min', stats.min);

    fprintf('%-12s | %20.4f\n', 'Max', stats.max);
    fprintf(fileID, '%-12s | %20.4f\n', 'Max', stats.max);
end

function plot_box_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val, 'z')
        plot_name = 'Z-Transformed Correlation Box Plot';
        y_label   = 'Fisher z-transform';
        suffix    = '_z';
    else
        plot_name = 'Correlation Box Plot';
        y_label   = 'Pearson Correlation (r)';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible', 'off', 'Position', [100,100,700,600]);
    h = gca;

    if ~isempty(data)
        boxplot(h, data', 'Notch', 'on', 'Labels', {'All'});
        margin = (stats.max - stats.min) * 0.1;
        if margin == 0 || isnan(margin), margin = 0.001; end
        ylim(h, [stats.min - margin, stats.max + margin]);
        title(h, plot_name);
        ylabel(h, y_label);
        grid(h, 'on');
    else
        title(h, 'All Correlations (No Data)');
        set(h, 'XTick', [], 'YTick', []);
    end

    fn = fullfile(results_dir, ['boxplot' suffix]);
    saveas(fig, [fn '.png']);
    saveas(fig, [fn '.svg']);
    savefig(fig, [fn '.fig']);
    fprintf('Saved boxplot%s.{png,svg,fig}\n', suffix);
    fprintf(fileID, 'Saved boxplot%s.{png,svg,fig}\n', suffix);
    close(fig);
end

function plot_hist_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val, 'z')
        plot_name = 'Z-Transformed Correlation Histogram';
        x_label   = 'Fisher z-transform';
        suffix    = '_z';
    else
        plot_name = 'Correlation Histogram';
        x_label   = 'Pearson Correlation (r)';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible', 'off', 'Position', [100,100,700,500]);
    h = gca;

    if ~isempty(data)
        histogram(h, data, 'BinMethod', 'auto');
        margin = (stats.max - stats.min) * 0.1;
        if margin == 0 || isnan(margin), margin = 0.001; end
        xlim(h, [stats.min - margin, stats.max + margin]);
    end

    title(h, [plot_name ' (Zoomed)']);
    xlabel(h, x_label);
    ylabel(h, 'Count');
    grid(h, 'on');

    fn = fullfile(results_dir, ['hist' suffix]);
    saveas(fig, [fn '.png']);
    saveas(fig, [fn '.svg']);
    savefig(fig, [fn '.fig']);
    fprintf('Saved hist%s.{png,svg,fig}\n', suffix);
    fprintf(fileID, 'Saved hist%s.{png,svg,fig}\n', suffix);
    close(fig);
end

function plot_z_score_distributions(~, z_stats, results_dir, fileID)
    z_mean = z_stats.mean;
    z_std  = z_stats.std;
    if z_std == 0 || isnan(z_std), z_std = 1e-6; end

    SE    = z_stats.sem;
    df    = max(z_stats.N - 1, 1);
    tcrit = tinv(0.975, df);
    CI    = [z_mean - tcrit*SE, z_mean + tcrit*SE];

    x = linspace(z_mean - 4*z_std, z_mean + 4*z_std, 1000);
    if ~isfinite(x(1)) || ~isfinite(x(end)) || numel(x) < 2 || x(1) == x(2)
        x = linspace(-3, 3, 1000);
    end

    pdf_data = normpdf(x, z_mean, z_std);

    fig = figure('Name','Z-Transformed Distributions', ...
                 'Visible','off', 'Position',[100,100,1000,600]);
    h = gca; hold(h, 'on');
    ph = plot(h, x, pdf_data, 'LineWidth', 2);
    yl = ylim(h);
    line(h, [z_mean z_mean], yl, 'LineStyle','--', 'LineWidth',1.5);
    line(h, [CI(1) CI(1)], yl, 'LineWidth',2);
    line(h, [CI(2) CI(2)], yl, 'LineWidth',2);
    ylim(h, yl);
    title(h,'Normal Distribution of Fisher z');
    xlabel(h,'z');
    ylabel(h,'PDF');
    legend(h, ph, {sprintf('PDF (mu=%.2f, std=%.2f)', z_mean, z_std)}, ...
           'Location','southeastoutside');
    grid(h,'on');
    hold(h,'off');

    fn = fullfile(results_dir,'z_score_distributions');
    saveas(fig, [fn '.png']);
    saveas(fig, [fn '.svg']);
    savefig(fig, [fn '.fig']);
    fprintf('Saved z_score_distributions.{png,svg,fig}\n');
    fprintf(fileID, 'Saved z_score_distributions.{png,svg,fig}\n');
    close(fig);
end

function print_raw_correlation_data_noID(all_results_table, all_correlations, fileID)
    fprintf('All correlations calculated.\n\n');
    fprintf(fileID, 'All correlations calculated.\n\n');

    fprintf('\nFull Experiment Summary\n');
    fprintf(fileID, '\nFull Experiment Summary\n');

    try
        unique_paths = unique(all_results_table.FilePath);
        fprintf('  Total Unique Traces: %d\n', numel(unique_paths));
        fprintf(fileID,'  Total Unique Traces: %d\n', numel(unique_paths));
        total_pairs = height(all_results_table);
        fprintf('  Total Pairwise Comparisons: %d\n', total_pairs);
        fprintf(fileID,'  Total Pairwise Comparisons: %d\n', total_pairs);
    catch
    end

    fprintf('\nAll Correlations (N=%d)\n', numel(all_correlations));
    fprintf(fileID, '\nAll Correlations (N=%d)\n', numel(all_correlations));

    for i = 1:numel(all_correlations)
        fprintf('%f  ', all_correlations(i));
        fprintf(fileID, '%f  ', all_correlations(i));
        if mod(i,5) == 0
            fprintf('\n');
            fprintf(fileID, '\n');
        end
    end
    fprintf('\n');
    fprintf(fileID, '\n');

    fprintf('\nFull Results Table\n');
    fprintf(fileID, '\nFull Results Table\n');
    try
        T_str = evalc('disp(all_results_table)');
        fprintf('%s\n', T_str);
        fprintf(fileID, '%s\n', T_str);
    catch
    end
end

function debug_dump_levels(tbl, ~, fileID)
    fprintf('\n[DEBUG] Levels present in tbl vs expected\n');
    fprintf(fileID, '\n[DEBUG] Levels present in tbl vs expected\n');

    sites = unique(tbl.Site);
    days  = unique(tbl.Day);
    times = unique(tbl.Time);

    fprintf('  Sites found: %s\n', strjoin(string(sites), ', '));
    fprintf('  Days  found: %s\n', strjoin(string(days),  ', '));
    fprintf('  Times found: %s\n', strjoin(string(times), ', '));

    fprintf(fileID, '  Sites found: %s\n', strjoin(string(sites), ', '));
    fprintf(fileID, '  Days  found: %s\n', strjoin(string(days),  ', '));
    fprintf(fileID, '  Times found: %s\n', strjoin(string(times), ', '));
end

% run_oneway_rm_day - 1-way repeated-measures ANOVA: Day factor (MON vs WED).
%
%   Averages correlations over Site within each (Time, Day) cell, then pivots
%   the table so Time is the repeated-measures subject and Day (MON, WED) is
%   the within-subject factor. The ANOVA is performed on Fisher z-transformed
%   correlations. Reports F, df, p, Greenhouse-Geisser corrected p, and partial
%   eta-squared. Generates and saves a dot-and-whisker plot of level means with
%   95% CI on the r-scale.
%
%   Inputs:
%     tbl    - Results table with Correlation, Site, Day, and Time columns
%     config - Configuration struct (provides allowed levels and results_dir)
%     fileID - Log file identifier for tee-printing
%
function run_oneway_rm_day(tbl, config, fileID)
    T = tbl;
    T.Site = upper(string(T.Site));
    T.Day  = upper(string(T.Day));
    T.Time = upper(string(T.Time));

    keep_days  = intersect(upper(string(config.allowed_days)),  unique(T.Day));
    keep_sites = intersect(upper(string(config.allowed_sites)), unique(T.Site));
    keep_times = intersect(upper(string(config.allowed_times)), unique(T.Time));

    T = T(ismember(T.Day, keep_days) & ...
          ismember(T.Site, keep_sites) & ...
          ismember(T.Time, keep_times), :);

    if height(T) == 0
        fprintf('  [DAY] No data after filtering.\n');
        fprintf(fileID, '  [DAY] No data after filtering.\n');
        return;
    end

    rz = max(min(T.Correlation,0.9999999), -0.9999999);
    T.z = atanh(rz);

    % Mean z per Time x Day (averaging over sites)
    G = groupsummary(T, {'Time','Day'}, 'mean', 'z');
    W = unstack(G, 'mean_z', 'Day', 'GroupingVariables', 'Time'); % MON, WED

    if ~ismember('MON', W.Properties.VariableNames), W.MON = NaN(height(W),1); end
    if ~ismember('WED', W.Properties.VariableNames), W.WED = NaN(height(W),1); end

    W = W(:, {'Time','MON','WED'});
    W = rmmissing(W, 'DataVariables', {'MON','WED'});

    nSubj = height(W);
    if nSubj < 2
        fprintf('  [DAY] Not enough subjects for RM-ANOVA.\n');
        fprintf(fileID, '  [DAY] Not enough subjects for RM-ANOVA.\n');
        debug_dump_levels(T, config, fileID);
        return;
    end

    if config.SCR_CFG_LIST_SUBJECT_KEYS
        fprintf('[DAY] Subjects (Time) used:\n');
        fprintf(fileID, '[DAY] Subjects (Time) used:\n');
        disp(W.Time);
        fprintf(fileID, '%s\n', evalc('disp(W.Time)'));
    end

    Y = W{:, {'MON','WED'}};
    tY = array2table(Y, 'VariableNames', {'MON','WED'});
    within = table(categorical({'MON';'WED'}), 'VariableNames', {'Day'});

    % Assumption checks
    try
        y  = Y(:);
        id = repelem((1:nSubj)', 2);
        mu_id = grpstats(y, id, 'mean');
        resid = y - mu_id(id);
        [~, pNorm] = adtest((resid-mean(resid))/std(resid));
        g = [ones(nSubj,1); 2*ones(nSubj,1)]; % MON=1, WED=2
        pLev = vartestn(y, g, 'TestType','LeveneAbsolute','Display','off');
        fprintf('  [DAY] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
        fprintf(fileID, '  [DAY] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
    catch
    end

    % 1-way RM-ANOVA
    rm  = fitrm(tY, 'MON-WED ~ 1', 'WithinDesign', within);
    ran = ranova(rm, 'WithinModel','Day');

    rn  = string(ran.Properties.RowNames);
    idx = find(rn=="(Intercept):Day",1);
    if isempty(idx), idx = find(contains(rn,'Day'),1); end
    errRow = find(contains(rn,'Error(Day)'),1);

    df1  = ran.DF(idx);
    df2  = ran.DF(errRow);
    Fval = ran.F(idx);
    p    = ran.pValue(idx);

    if ismember('pValueGG', ran.Properties.VariableNames) && ~isnan(ran.pValueGG(idx))
        pGG = ran.pValueGG(idx);
    else
        pGG = p;
    end

    SS_eff = ran.SumSq(idx);
    SS_err = ran.SumSq(errRow);
    eta_p2 = SS_eff / (SS_eff + SS_err);

    % Print compact ANOVA summary table (ASCII)
    fprintf('    RM-ANOVA Summary (Day)\n');
    fprintf('    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf('    Day      | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    fprintf(fileID, '    RM-ANOVA Summary (Day)\n');
    fprintf(fileID, '    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf(fileID, '    Day      | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    % Level means + CI on r-scale and dot-whisker plot
    levels = {'MON','WED'};
    [mu_r, lo_r, hi_r] = print_level_means_ci_from_Y('Day', levels, Y, fileID);
    plot_dotwhisker_r(mu_r, lo_r, hi_r, levels, config.results_dir, 'day', ...
                      'Estimates with 95% CI - DAY (MON vs WED)');
end

% run_oneway_rm_site - 1-way repeated-measures ANOVA: Site factor (SBOX vs SPSU).
%
%   Averages correlations over Day within each (Time, Site) cell, then pivots
%   the table so Time is the repeated-measures subject and Site (SBOX, SPSU)
%   is the within-subject factor. The ANOVA is performed on Fisher z-transformed
%   correlations. Reports F, df, p, Greenhouse-Geisser corrected p, and partial
%   eta-squared. Generates and saves a dot-and-whisker plot of level means with
%   95% CI on the r-scale.
%
%   Inputs:
%     tbl    - Results table with Correlation, Site, Day, and Time columns
%     config - Configuration struct (provides allowed levels and results_dir)
%     fileID - Log file identifier for tee-printing
%
function run_oneway_rm_site(tbl, config, fileID)
    T = tbl;
    T.Site = upper(string(T.Site));
    T.Day  = upper(string(T.Day));
    T.Time = upper(string(T.Time));

    keep_days  = intersect(upper(string(config.allowed_days)),  unique(T.Day));
    keep_sites = intersect(upper(string(config.allowed_sites)), unique(T.Site));
    keep_times = intersect(upper(string(config.allowed_times)), unique(T.Time));

    T = T(ismember(T.Day, keep_days) & ...
          ismember(T.Site, keep_sites) & ...
          ismember(T.Time, keep_times), :);

    if height(T) == 0
        fprintf('  [SITE] No data after filtering.\n');
        fprintf(fileID, '  [SITE] No data after filtering.\n');
        return;
    end

    rz = max(min(T.Correlation,0.9999999), -0.9999999);
    T.z = atanh(rz);

    % Mean z per Time x Site (averaging over days)
    G = groupsummary(T, {'Time','Site'}, 'mean', 'z');
    W = unstack(G, 'mean_z', 'Site', 'GroupingVariables', 'Time'); % SBOX, SPSU

    if ~ismember('SBOX', W.Properties.VariableNames), W.SBOX = NaN(height(W),1); end
    if ~ismember('SPSU', W.Properties.VariableNames), W.SPSU = NaN(height(W),1); end

    W = W(:, {'Time','SBOX','SPSU'});
    W = rmmissing(W, 'DataVariables', {'SBOX','SPSU'});

    nSubj = height(W);
    if nSubj < 2
        fprintf('  [SITE] Not enough subjects for RM-ANOVA.\n');
        fprintf(fileID, '  [SITE] Not enough subjects for RM-ANOVA.\n');
        debug_dump_levels(T, config, fileID);
        return;
    end

    if config.SCR_CFG_LIST_SUBJECT_KEYS
        fprintf('[SITE] Subjects (Time) used:\n');
        fprintf(fileID, '[SITE] Subjects (Time) used:\n');
        disp(W.Time);
        fprintf(fileID, '%s\n', evalc('disp(W.Time)'));
    end

    Y = W{:, {'SBOX','SPSU'}};
    tY = array2table(Y, 'VariableNames', {'SBOX','SPSU'});
    within = table(categorical({'SBOX';'SPSU'}), 'VariableNames', {'Site'});

    % Assumption checks
    try
        y  = Y(:);
        id = repelem((1:nSubj)', 2);
        mu_id = grpstats(y, id, 'mean');
        resid = y - mu_id(id);
        [~, pNorm] = adtest((resid-mean(resid))/std(resid));
        g = [ones(nSubj,1); 2*ones(nSubj,1)]; % SBOX=1, SPSU=2
        pLev = vartestn(y, g, 'TestType','LeveneAbsolute','Display','off');
        fprintf('  [SITE] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
        fprintf(fileID, '  [SITE] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
    catch
    end

    % 1-way RM-ANOVA
    rm  = fitrm(tY, 'SBOX-SPSU ~ 1', 'WithinDesign', within);
    ran = ranova(rm, 'WithinModel', 'Site');

    rn  = string(ran.Properties.RowNames);
    idx = find(rn=="(Intercept):Site",1);
    if isempty(idx), idx = find(contains(rn,'Site'),1); end
    errRow = find(contains(rn,'Error(Site)'),1);

    df1  = ran.DF(idx);
    df2  = ran.DF(errRow);
    Fval = ran.F(idx);
    p    = ran.pValue(idx);

    if ismember('pValueGG', ran.Properties.VariableNames) && ~isnan(ran.pValueGG(idx))
        pGG = ran.pValueGG(idx);
    else
        pGG = p;
    end

    SS_eff = ran.SumSq(idx);
    SS_err = ran.SumSq(errRow);
    eta_p2 = SS_eff / (SS_eff + SS_err);

    % Print compact ANOVA summary table (ASCII)
    fprintf('    RM-ANOVA Summary (Site)\n');
    fprintf('    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf('    Site     | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    fprintf(fileID, '    RM-ANOVA Summary (Site)\n');
    fprintf(fileID, '    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf(fileID, '    Site     | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    % Level means + CI on r-scale and dot-whisker plot
    levels = {'SBOX','SPSU'};
    [mu_r, lo_r, hi_r] = print_level_means_ci_from_Y('Site', levels, Y, fileID);
    plot_dotwhisker_r(mu_r, lo_r, hi_r, levels, config.results_dir, 'site', ...
                      'Estimates with 95% CI - SITE (SBOX vs SPSU)');
end

% print_level_means_ci_from_Y - Compute and print per-level means with 95% CI (r-scale).
%
%   Operates on the z-scale subject-by-level matrix Y. Computes per-level
%   mean and SE on the z-scale, derives the 95% CI using the t-critical value,
%   then back-transforms mean and CI bounds to r-scale via tanh. Also prints
%   both z-scale and r-scale summaries to the console and log.
%
%   Inputs:
%     tag    - Factor name string used in printed output (e.g., 'DAY')
%     names  - 1x2 cell array of level name strings (e.g., {'MON','WED'})
%     Y      - nSubj x 2 matrix of Fisher z values (one column per level)
%     fileID - Log file identifier for tee-printing
%
%   Outputs:
%     mu_r  - 1x2 vector of mean Pearson r values (one per level)
%     lo_r  - 1x2 vector of 95% CI lower bounds (r-scale)
%     hi_r  - 1x2 vector of 95% CI upper bounds (r-scale)
%
function [mu_r, lo_r, hi_r] = print_level_means_ci_from_Y(tag, names, Y, fileID)
    % Y: nSubj x 2, z-scale values for the two levels

    n    = size(Y,1);
    mu_z = mean(Y,1);
    sd_z = std(Y,0,1);
    se_z = sd_z / max(sqrt(n),1);
    df   = max(n-1,1);
    tcrit = tinv(0.975, df);

    lo_z = mu_z - tcrit*se_z;
    hi_z = mu_z + tcrit*se_z;

    mu_r = tanh(mu_z);
    lo_r = tanh(lo_z);
    hi_r = tanh(hi_z);

    fprintf('  [%s] Level means (z''): %s=%.4f [%.4f, %.4f],  %s=%.4f [%.4f, %.4f]\n', ...
        tag, names{1}, mu_z(1), lo_z(1), hi_z(1), ...
             names{2}, mu_z(2), lo_z(2), hi_z(2));

    fprintf('  [%s] Level 95%% CI (r): %s=[%.4f,%.4f],  %s=[%.4f,%.4f]\n', ...
        tag, names{1}, lo_r(1), hi_r(1), ...
             names{2}, lo_r(2), hi_r(2));

    if fileID > 1
        fprintf(fileID, '  [%s] Level means (z''): %s=%.4f [%.4f, %.4f],  %s=%.4f [%.4f, %.4f]\n', ...
            tag, names{1}, mu_z(1), lo_z(1), hi_z(1), ...
                 names{2}, mu_z(2), lo_z(2), hi_z(2));
        fprintf(fileID, '  [%s] Level 95%% CI (r): %s=[%.4f,%.4f],  %s=[%.4f,%.4f]\n', ...
            tag, names{1}, lo_r(1), hi_r(1), ...
                 names{2}, lo_r(2), hi_r(2));
    end
end

function plot_dotwhisker_r(mu_r, ci_lo_r, ci_hi_r, levels, results_dir, tag, title_str)
    % Dot-and-whisker plot for means 95% CIs on r-scale.

    if nargin < 7 || isempty(title_str)
        title_str = ['Estimates with 95% CI - ' upper(tag)];
    end

    mu_r    = mu_r(:)';
    ci_lo_r = ci_lo_r(:)';
    ci_hi_r = ci_hi_r(:)';
    K       = numel(mu_r);
    x       = 1:K;

    % Simple color scheme
    mainColor    = [63  81 181] / 255;
    dotFaceColor = [26  35 126] / 255;
    dotEdgeColor = [1 1 1];

    fig = figure('Name',['DotWhisker_' tag], ...
                 'Visible','off', ...
                 'Position',[120,120,820,560]);
    ax = gca;
    hold(ax,'on');

    for i = 1:K
        line(ax, [x(i) x(i)], [ci_lo_r(i) ci_hi_r(i)], ...
            'Color', mainColor, ...
            'LineWidth', 2.0);
        plot(ax, x(i), mu_r(i), 'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', dotFaceColor, ...
            'MarkerEdgeColor', dotEdgeColor, ...
            'LineWidth', 1.2);
    end

    set(ax, 'XLim', [0.5 K+0.5], ...
            'XTick', x, ...
            'XTickLabel', levels, ...
            'FontName', 'Times', ...
            'FontSize', 10, ...
            'Box', 'off', ...
            'Layer', 'top', ...
            'YGrid', 'on', ...
            'LineWidth', 0.8);

    ylabel(ax, 'Pearson r', 'FontSize', 10);

    vals = [ci_lo_r(:); ci_hi_r(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vals = [0.95 1.0];
    end
    ymin = max(0, min(vals) - 0.001);
    ymax = min(1.0, max(vals) + 0.001);
    if ymin >= ymax
        ymin = max(0, ymin - 0.01);
        ymax = min(1.0, ymax + 0.01);
    end
    ylim(ax, [ymin ymax]);

    step = 0.001;
    yticks = ceil(ymin/step)*step : step : floor(ymax/step)*step;
    if numel(yticks) < 2
        yticks = [ymin ymax];
    end
    set(ax, 'YTick', yticks);

    title(ax, title_str, 'FontWeight','normal', 'FontSize',10);
    grid(ax,'on');
    hold(ax,'off');

    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end

    fn = fullfile(results_dir, ['dotwhisker_' tag]);
    saveas(fig, [fn '.png']);
    saveas(fig, [fn '.svg']);
    savefig(fig, [fn '.fig']);
    close(fig);

    fprintf('    Saved dotwhisker_%s.{png,svg,fig}\n', tag);
end
