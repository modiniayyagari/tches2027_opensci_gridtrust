%% ENF Temporal Reliability and Replay-Resilience Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Quantifies how stable the ENF-based correlation between a DC-powered FPGA
%   and the AC mains reference remains across time of day, day of week, and week
%   at a single grid location (TREND dataset, US 60 Hz grid). 
%
%   Additionally runs an older-reference stale-replay attack to evaluate threshold-based rejection
%   of replayed ambient traces against earlier verifier reference windows.
%
% Analysis Steps:
%   1. Scan exp_inputs/TREND/ for matched mains_pow_trace_ac.wav and
%      fpga_em_trace_dc.wav pairs. Parse week, day-type, day, time-of-day,
%      and leaf-folder metadata from the path tokens.
%   2. Compute Pearson correlations (r) using proc_enf_analysis.
%   3. Apply the Fisher z = atanh(r) transform for variance stabilization.
%   4. Report descriptive statistics (N, mean, SD, 95% CI) for r and z.
%   5. Run within-condition repeated-measures ANOVA tests:
%        - Time-of-day: EMRN vs MORN vs AFTN vs EVEN (subjects = Week x Day)
%        - Day-of-week: WED vs THU vs SAT vs SUN (subjects = Week x TimeOfDay)
%        - Week-to-week: WK01 vs WK02 paired t-test (subjects = Day x TimeOfDay)
%   6. Apply Holm-Bonferroni correction across the two RM tests.
%   7. Generate a combined three-panel dot-and-whisker figure for all factors.
%   8. Run an older-reference stale-replay attack: compare each target ambient
%      trace against all strictly older mains-reference windows and report
%      pooled authentic-vs-replay correlation statistics.
%
% Inputs:
%   exp_inputs/TREND/ - Artifact dataset with WK01/ and WK02/ sub-trees, each
%                       containing WDAY/ and WEND/ sub-trees with T01-T05 leaf
%                       folders holding the paired .wav trace files.
%
% Outputs:
%   exp_results/enf_temporal_reliability/ - All plot files (.png, .svg, .fig)
%   exp_results/replay_attack_resilience/ - Replay attack CSV and ENF overlay
%   exp_logs/<script_name>_log.txt        - Full text log of all output
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

%% 1) Configuration
% Set analysis flags, STFT parameters, ENF extraction settings, and
% replay attack target specifications for the TREND dataset.
config = struct();

% Resolve the artifact data root relative to this script location.
script_dir    = fileparts(mfilename('fullpath'));
artifact_root = fileparts(script_dir);
config.baseDir = fullfile(artifact_root, 'exp_inputs', 'TREND');

% WAV base names used in each TREND trial folder.
config.file_1_name_base = 'mains_pow_trace_ac'; % mains reference
config.file_2_name_base = 'fpga_em_trace_dc';   % sensed ambient trace

% Analysis and reporting flags.
config.SHOW_RAW_CORRELATIONS      = false;   % print each r value and table row
config.SHOW_DESCRIPTIVE_STATS     = true;    % write descriptive stats and plots
config.RUN_INFERENTIAL_TESTS      = true;    % run TimeOfDay, Day, and Week tests
config.APPLY_HOLM_CORRECTION      = true;    % correct TimeOfDay and Day p-values
config.LIST_RM_SUBJECTS           = false;   % list subjects retained in RM tests
config.SAVE_QQ_PLOTS              = false;   % save residual/difference Q-Q plots
config.PLOT_COMBINED_INFERENTIAL  = true;    % export one three-panel summary plot

% Replay attack target traces.
% Each row is {Week, DayType, Day, TimeOfDay, LeafFolder}. THU is normalized
% to the dataset's THUR folder name when needed.
config.replay_target_specs = {
    'WK01','WDAY','THU', 'AFTN','T01';
    'WK01','WEND','SAT', 'EMRN','T02';
    'WK01','WEND','SAT', 'EVEN','T03';
    'WK01','WEND','SUN', 'MORN','T04';
    'WK01','WEND','SUN', 'EVEN','T05';
    'WK02','WDAY','WED', 'EMRN','T01';
    'WK02','WEND','SAT', 'MORN','T02';
    'WK02','WEND','SUN', 'AFTN','T03';
};
config.replay_overlay_target_spec = {'WK02','WEND','SUN','AFTN','T02'};

% --- Leaf folder filter (limit to specific Txx folders) ---
config.FILTER_BY_LEAF_FOLDER      = true;
config.allowed_leaf_folders       = {'T01','T02','T03','T04','T05'};

% --- Spectrogram / ENF extraction parameters ---
config.nominal_freq_arr = [50 60];      % [50Hz, 60Hz]
frame_size_arr          = (1:12)*1000;
config.frame_size       = frame_size_arr(8);  % 8000 ms
nfft_arr                = 2.^(10:20);
config.nfft             = nfft_arr(6);        % 2^15
overlap_size_arr        = 0:0.1:0.9;
config.overlap_size     = overlap_size_arr(1)*config.frame_size;

% Trace-specific parameters (for proc_enf_analysis)
nominal_freq = config.nominal_freq_arr(2); % 60 Hz
harmonics = (1:7) * nominal_freq;

config.harmonics_arr_1                        = harmonics;
config.trace_1_freq_est_method                = 1;
config.trace_1_est_freq                       = nominal_freq;
config.trace_1_freq_est_spec_comb_harmonics   = harmonics;
config.trace_1_plot_title                     = config.file_1_name_base;

config.harmonics_arr_2                        = harmonics;
config.trace_2_freq_est_method                = 1;
config.trace_2_est_freq                       = nominal_freq;
config.trace_2_freq_est_spec_comb_harmonics   = harmonics;
config.trace_2_plot_title                     = config.file_2_name_base;

%% 2) Setup: Logging and Colormap
% Open the log file, create output directories, and apply the custom
% diverging colormap used by all saved figures in this script.
[fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders();
config.results_dir = results_dir;
setup_colormap();

fprintf('Starting ENF Temporal Reliability Analysis...\n');
fprintf(fileID, 'Starting ENF Temporal Reliability Analysis. Log: %s\n\n', log_filename);

%% 3) Main Analysis: Batch ENF Correlation
% Scan all TREND trace pairs, call proc_enf_analysis on each, and
% collect the Pearson correlation results into a labelled results table.
results_tbl      = build_results_table(config, fileID);
all_correlations = results_tbl.Correlation;


fprintf('Analysis complete. N=%d correlations.\n\n', numel(all_correlations));
fprintf(fileID, 'Analysis complete. N=%d correlations.\n\n', numel(all_correlations));

%% 4) Raw Results (Optional)
% Optionally print every individual correlation value and the full
% results table to the console and log file.
if config.SHOW_RAW_CORRELATIONS
    print_raw_correlations(results_tbl, all_correlations, fileID);
end

%% 5) Descriptive Statistics and Plots
% Compute mean, SD, 95% CI, and Fisher z statistics over all collected
% correlations; save box plots, histograms, and z-distribution figures.
[r_stats, z_data, z_stats] = calculate_basic_stats(all_correlations);

if config.SHOW_DESCRIPTIVE_STATS
    print_descriptive_stats_results(r_stats, z_stats, fileID);
    try
        plot_box_simple(all_correlations, r_stats, results_dir, fileID, 'r');
        plot_hist_simple(all_correlations, r_stats, results_dir, fileID, 'r');
        plot_box_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_hist_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_z_score_distributions(z_data, z_stats, results_dir, fileID);
    catch MEplots
        fprintf('Plot error: %s\n', MEplots.message);
        fprintf(fileID, 'Plot error: %s\n', MEplots.message);
    end
else
    fprintf('\nSTEP 1 (Descriptive) skipped by config.\n');
    fprintf(fileID, '\nSTEP 1 (Descriptive) skipped by config.\n');
end

%% 6) Inferential Statistics: Time-of-Day, Day, and Week
% Run within-condition RM-ANOVA tests for time-of-day and day-of-week effects,
% a paired t-test for week-to-week stability, Holm-Bonferroni correction across
% the two RM tests, and export a combined three-panel dot-whisker summary figure.
if config.RUN_INFERENTIAL_TESTS
    fprintf('\n\nSTEP 2: Inferential Tests (within-condition repeated measures)\n');
    fprintf(fileID, '\n\nSTEP 2: Inferential Tests (within-condition repeated measures)\n');

    try
        week_tbl = results_tbl(results_tbl.Week ~= "NA", :);

        % Time-of-day RM-ANOVA (z-scale; report on r-scale)
        out_time = run_rm_anova(week_tbl, ...
            {'Week','Day'}, 'TimeOfDay', ...
            {'EMRN','MORN','AFTN','EVEN'}, fileID, config);

        % Day-of-week RM-ANOVA
        out_day = run_rm_anova(week_tbl, ...
            {'Week','TimeOfDay'}, 'Day', ...
            {'WED','THU','SAT','SUN'}, fileID, config);

        % Week-to-week paired comparison
        out_week = run_week_paired(week_tbl, fileID, config);

        % Holm-Bonferroni across the two RM tests (GG p-values)
        if config.APPLY_HOLM_CORRECTION
            apply_holm_two_tests(out_time, out_day, fileID);
        end

        % Combined figure with all three dot-whisker panels
        if config.PLOT_COMBINED_INFERENTIAL
            try
                plot_combined_inferential_boxplots(out_time, out_day, out_week, results_dir);
            catch MEc
                fprintf('  [WARN] Combined inferential plot failed: %s\n', MEc.message);
                fprintf(fileID, '  [WARN] Combined inferential plot failed: %s\n', MEc.message);
            end
        end

    catch MEopt
        fprintf('  [ERROR] Inferential tests failed: %s\n', MEopt.message);
        fprintf(fileID, '  [ERROR] Inferential tests failed: %s\n', MEopt.message);
    end
else
    fprintf('\n\nSTEP 2 skipped by config.\n');
    fprintf(fileID, '\n\nSTEP 2 skipped by config.\n');
end

%% 7) Replay Attack Resilience Analysis (Older-Reference Stale-Replay Test)
% Evaluate how well a fixed decision threshold rejects a stale ambient trace
% when it is compared against strictly older verifier reference windows,
% modelling an older-reference stale-replay attack (captured at time T replayed
% against verifier windows at time T'' < T).
fprintf('\n\nSTEP 3: Replay Attack Resilience (Older Reference Stale Replay)...\n');
fprintf(fileID, '\n\nSTEP 3: Replay Attack Resilience (Older Reference Stale Replay)...\n');
try
    run_replay_attack_analysis(config, fileID);
catch ME_replay
    fprintf('  [ERROR] Replay attack analysis failed: %s\n', ME_replay.message);
    fprintf(fileID, '  [ERROR] Replay attack analysis failed: %s\n', ME_replay.message);
end

fprintf('\n\nAll analysis complete. Log saved.\n');
fprintf(fileID, '\n\nAll analysis complete.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup_logging_and_folders - Create output directories and open the log file.
%
%   Opens (or overwrites) a timestamped log file under exp_logs/ and creates
%   the exp_results/enf_temporal_reliability/ results directory. Returns a
%   cleanup object that closes the file descriptor when it goes out of scope.
%
%   Inputs:
%     (none)
%
%   Outputs:
%     fileID     - File descriptor for the open log file (from fopen).
%     results_dir - Path string of the results output directory.
%     log_filename - Full path of the log file that was opened.
%     cleanupObj  - onCleanup object that closes fileID on function exit.
%
function [fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders()
    log_dir     = 'exp_logs';
    results_dir = 'exp_results/enf_temporal_reliability';
    if ~exist(log_dir, 'dir'), mkdir(log_dir); end
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    [~, script_name] = fileparts(mfilename('fullpath'));
    if isempty(script_name)
        script_name = 'enf_temporal_reliability';
    end

    log_filename = fullfile(log_dir, [script_name, '_log.txt']);
    if exist(log_filename, 'file'), delete(log_filename); end

    fileID = fopen(log_filename, 'w');
    if fileID == -1, error('Could not open %s', log_filename); end

    cleanupObj = onCleanup(@() fclose(fileID));
end

function setup_colormap()
    numColors = 256;
    pos_dark  = [144, 202, 249] / 255;
    neu_color = [255, 255, 255] / 255;
    neg_dark  = [229, 115, 115] / 255;
    cmap_negative = [linspace(neg_dark(1), neu_color(1), numColors/2)', ...
                     linspace(neg_dark(2), neu_color(2), numColors/2)', ...
                     linspace(neg_dark(3), neu_color(3), numColors/2)'];
    cmap_positive = [linspace(neu_color(1), pos_dark(1), numColors/2)', ...
                     linspace(neu_color(2), pos_dark(2), numColors/2)', ...
                     linspace(neu_color(3), pos_dark(3), numColors/2)'];
    set(0,'DefaultFigureColormap',[cmap_negative; cmap_positive]);
end

% build_results_table - Scan the TREND dataset and compute Pearson correlations.
%
%   Recursively finds every fpga_em_trace_dc.wav file under config.baseDir,
%   matches each with the co-located mains_pow_trace_ac.wav reference, calls
%   proc_enf_analysis, and records the correlation together with parsed path
%   metadata (DayType, TimeOfDay, Week, Day, LeafFolder).
%
%   Inputs:
%     config  - Configuration struct with fields: baseDir, file_1_name_base,
%               file_2_name_base, FILTER_BY_LEAF_FOLDER, allowed_leaf_folders,
%               nfft, frame_size, overlap_size, harmonics_arr_1/2,
%               nominal_freq_arr, trace_1/2_freq_est_method, _est_freq,
%               _freq_est_spec_comb_harmonics, _plot_title.
%     fileID  - Log file descriptor (from setup_logging_and_folders).
%
%   Outputs:
%     tbl - MATLAB table with columns: Correlation, FilePath, DayType,
%           TimeOfDay, Week, Day, LeafFolder.
%
function tbl = build_results_table(config, fileID)
    % Scan for sensed files, match with mains, run ENF correlation.
    varTypes = {'double','string','string','string','string','string','string'};
    varNames = {'Correlation','FilePath','DayType','TimeOfDay','Week','Day','LeafFolder'};

    fprintf('INFO: Searching for trace pairs in: %s\n', config.baseDir);
    fprintf(fileID, 'INFO: Searching for trace pairs in: %s\n', config.baseDir);

    sensed_files = dir(fullfile(config.baseDir, '**', [config.file_2_name_base, '.wav']));
    fprintf('INFO: Found %d sensed files.\n', numel(sensed_files));
    fprintf(fileID, 'INFO: Found %d sensed files.\n', numel(sensed_files));

    tbl = table('Size',[numel(sensed_files), numel(varNames)], ...
                'VariableTypes',varTypes, 'VariableNames',varNames);
    rowCount = 0;

    hasToken = @(p, tok) ~isempty(regexp(p, ['(^|[\\/_.-])' tok '([\\/_.-]|$)'], 'once'));

    for k = 1:numel(sensed_files)
        currentFolder = sensed_files(k).folder;

        [~, leaf] = fileparts(currentFolder);

        % Optional leaf folder filter
        if config.FILTER_BY_LEAF_FOLDER
            if ~ismember(leaf, config.allowed_leaf_folders)
                fprintf('\n--- Skip %d/%d (Leaf: %s not in allowed list) ---\n', ...
                        k, numel(sensed_files), leaf);
                fprintf(fileID, '\n--- Skip %d/%d (Leaf: %s not in allowed list) ---\n', ...
                        k, numel(sensed_files), leaf);
                continue;
            end
        end

        f2 = fullfile(currentFolder, [config.file_2_name_base, '.wav']); % sensed
        f1 = fullfile(currentFolder, [config.file_1_name_base, '.wav']); % reference

        if ~isfile(f1)
            fprintf('WARN: Missing reference WAV at %s\n', currentFolder);
            fprintf(fileID, 'WARN: Missing reference WAV at %s\n', currentFolder);
            continue;
        end

        % Parse labels from folder path
        dayType = 'NA';
        if hasToken(currentFolder, 'WDAY')
            dayType = 'WDAY';
        elseif hasToken(currentFolder, 'WEND')
            dayType = 'WEND';
        end

        day = 'NA';
        if hasToken(currentFolder, 'MON')
            day = 'MON';
        elseif hasToken(currentFolder, 'TUE')
            day = 'TUE';
        elseif hasToken(currentFolder, 'WED')
            day = 'WED';
        elseif hasToken(currentFolder,'THU') || hasToken(currentFolder,'THUR')
            day = 'THU';
        elseif hasToken(currentFolder, 'FRI')
            day = 'FRI';
        elseif hasToken(currentFolder, 'SAT')
            day = 'SAT';
        elseif hasToken(currentFolder, 'SUN')
            day = 'SUN';
        end

        timeOfDay = 'NA';
        if hasToken(currentFolder, 'EMRN')
            timeOfDay = 'EMRN';
        elseif hasToken(currentFolder, 'MORN')
            timeOfDay = 'MORN';
        elseif hasToken(currentFolder, 'AFTN')
            timeOfDay = 'AFTN';
        elseif hasToken(currentFolder, 'EVEN')
            timeOfDay = 'EVEN';
        end

        week = 'NA';
        if hasToken(currentFolder, 'WK01')
            week = 'WK01';
        elseif hasToken(currentFolder, 'WK02')
            week = 'WK02';
        end

        try
            r_corr = proc_enf_analysis( ...
                    f1, f2, ...
                    config.nfft, config.frame_size, config.overlap_size, ...
                    config.harmonics_arr_1, config.nominal_freq_arr(2), ...
                    config.harmonics_arr_2, config.nominal_freq_arr(2), ...
                    config.trace_1_freq_est_method, config.trace_1_est_freq, config.trace_1_freq_est_spec_comb_harmonics, ...
                    config.trace_2_freq_est_method, config.trace_2_est_freq, config.trace_2_freq_est_spec_comb_harmonics, ...
                    config.trace_1_plot_title, config.trace_2_plot_title, false);

            rowCount = rowCount + 1;
            tbl(rowCount,:) = {r_corr, string(f2), string(dayType), string(timeOfDay), string(week), string(day), string(leaf)};

            fprintf('\nINFO: %d/%d corr=%.6f', ...
                k, numel(sensed_files), r_corr);
            fprintf(fileID, '\nINFO: %d/%d corr=%.6f', ...
                k, numel(sensed_files), r_corr);
        catch ME
            fprintf('ERROR: Analysis failed for pair %d: %s\n', k, ME.message);
            fprintf(fileID, 'ERROR: Analysis failed for pair %d: %s\n', k, ME.message);
        end
    end

    tbl = tbl(1:rowCount, :);
    fprintf('\nINFO: Batch complete. N=%d\n', height(tbl));
    fprintf(fileID, '\nINFO: Batch complete. N=%d\n', height(tbl));
end


% calculate_basic_stats - Compute descriptive statistics on r and Fisher z.
%
%   Clips raw Pearson r values to (-1, 1), applies the Fisher z = atanh(r)
%   transform, and returns descriptive statistics structs for both scales.
%
%   Inputs:
%     all_correlations - Vector of raw Pearson r values.
%
%   Outputs:
%     r_stats - Descriptive stats struct on the r scale (see sub_calc_stats).
%     z_data  - Fisher z-transformed values (same length as all_correlations).
%     z_stats - Descriptive stats struct on the z scale (see sub_calc_stats).
%
function [r_stats, z_data, z_stats] = calculate_basic_stats(all_correlations)
    r_stats = sub_calc_stats(all_correlations);

    r_safe = all_correlations;
    r_safe(r_safe >= 1.0)  = 0.9999999;
    r_safe(r_safe <= -1.0) = -0.9999999;

    z_data  = atanh(r_safe);
    z_stats = sub_calc_stats(z_data);
end

% sub_calc_stats - Compute N, mean, median, SD, SEM, CI, skewness, kurtosis.
%
%   Strips non-finite values and returns a scalar struct with common
%   descriptive statistics. Returns NaN fields for empty or all-NaN inputs.
%
%   Inputs:
%     x - Numeric vector (any orientation).
%
%   Outputs:
%     stats - Struct with fields: N, mean, median, std, min, max, var,
%             skew, kurt, sem, ci95_margin.
%
function stats = sub_calc_stats(x)
    x = x(isfinite(x));
    if isempty(x)
        stats = struct('N',0,'mean',NaN,'median',NaN,'std',NaN,'min',NaN,'max',NaN, ...
                       'var',NaN,'skew',NaN,'kurt',NaN,'sem',NaN,'ci95_margin',NaN);
        return;
    end
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
    df           = max(stats.N - 1, 1);
    tcrit        = tinv(0.975, df);
    stats.ci95_margin = tcrit * stats.sem;
end

% print_descriptive_stats_results - Print r-scale and Fisher z descriptive stats.
%
%   Reports Step 1A (raw Pearson r statistics) and Step 1B (Fisher z statistics
%   with back-transformed mean and 95% CI) to both the console and the log file.
%
%   Inputs:
%     r_stats - Stats struct for raw Pearson r (from sub_calc_stats).
%     z_stats - Stats struct for Fisher z values (from sub_calc_stats).
%     fileID  - Log file descriptor.
%
%   Outputs:
%     (none) - Prints formatted tables to console and log file.
%
function print_descriptive_stats_results(r_stats, z_stats, fileID)
    % STEP 1A: r-scale
    fprintf('\nSTEP 1A: Descriptive Results (r)\n');
    fprintf(fileID, '\nSTEP 1A: Descriptive Results (r)\n');
    print_stats_table(r_stats, fileID, 'All Correlations (r)');

    % STEP 1B: Fisher z
    fprintf('\nSTEP 1B: Descriptive Results (z)\n');
    fprintf(fileID, '\nSTEP 1B: Descriptive Results (z)\n');

    % back-transformed mean and CI on r-scale
    z_mean        = z_stats.mean;
    mean_r_from_z = tanh(z_mean);
    ci_z          = [z_stats.mean - z_stats.ci95_margin, z_stats.mean + z_stats.ci95_margin];
    ci_r_from_z   = tanh(ci_z);

    print_stats_table(z_stats, fileID, 'All Correlations (z'')');
    fprintf('Mean (r from z-mean)   | %20.4f\n', mean_r_from_z);
    fprintf('95%% CI (r from z)      | [%0.4f, %0.4f]\n', ci_r_from_z(1), ci_r_from_z(2));
    fprintf(fileID, 'Mean (r from z-mean)   | %20.4f\n', mean_r_from_z);
    fprintf(fileID, '95%% CI (r from z)      | [%0.4f, %0.4f]\n', ci_r_from_z(1), ci_r_from_z(2));

    fprintf('STEP 1 complete.\n');
    fprintf(fileID, 'STEP 1 complete.\n');
end

function print_stats_table(stats, fileID, label)
    fprintf('%-20s | %20s\n', 'Statistic', label);
    fprintf(fileID, '%-20s | %20s\n', 'Statistic', label);

    fprintf('%-20s | %20d\n', 'Count (N)', stats.N);
    fprintf('%-20s | %20.4f\n', 'Mean', stats.mean);
    fprintf('%-20s | %20.4f\n', 'Median', stats.median);
    fprintf('%-20s | %20.4f\n', 'Std. Dev.', stats.std);
    fprintf('%-20s | %20.4f\n', 'Variance', stats.var);
    fprintf('%-20s | %20.4f\n', 'Skewness', stats.skew);
    fprintf('%-20s | %20.4f\n', 'Kurtosis', stats.kurt);
    fprintf('%-20s | %20.4f\n', 'Std. Error', stats.sem);
    fprintf('%-20s | [%0.4f, %0.4f]\n', '95% CI', ...
            stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin);
    fprintf('%-20s | %20.4f\n', 'Min', stats.min);
    fprintf('%-20s | %20.4f\n', 'Max', stats.max);

    fprintf(fileID, '%-20s | %20d\n', 'Count (N)', stats.N);
    fprintf(fileID, '%-20s | %20.4f\n', 'Mean', stats.mean);
    fprintf(fileID, '%-20s | %20.4f\n', 'Median', stats.median);
    fprintf(fileID, '%-20s | %20.4f\n', 'Std. Dev.', stats.std);
    fprintf(fileID, '%-20s | %20.4f\n', 'Variance', stats.var);
    fprintf(fileID, '%-20s | %20.4f\n', 'Skewness', stats.skew);
    fprintf(fileID, '%-20s | %20.4f\n', 'Kurtosis', stats.kurt);
    fprintf(fileID, '%-20s | %20.4f\n', 'Std. Error', stats.sem);
    fprintf(fileID, '%-20s | [%0.4f, %0.4f]\n', '95% CI', ...
            stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin);
    fprintf(fileID, '%-20s | %20.4f\n', 'Min', stats.min);
    fprintf(fileID, '%-20s | %20.4f\n', 'Max', stats.max);
end

% run_rm_anova - One-way repeated-measures ANOVA with Greenhouse-Geisser correction.
%
%   Aggregates Fisher z per subject (defined by subjectGroupingVars) and factor
%   level, fits a repeated-measures model via fitrm/ranova, reports F, p, p(GG),
%   and partial eta-squared, and saves a dot-whisker plot of level means on the
%   r scale with 95% confidence intervals.
%
%   Inputs:
%     tbl                - Results table with columns: Correlation, Week, Day,
%                          TimeOfDay, plus any grouping columns.
%     subjectGroupingVars - Cell array of column names defining RM subjects
%                           (e.g., {'Week','Day'} for TimeOfDay factor).
%     factorName         - Name of the within-subjects factor column (string).
%     levels             - Cell array of level labels in analysis order.
%     fileID             - Log file descriptor.
%     config             - Config struct (needs results_dir, SAVE_QQ_PLOTS,
%                          LIST_RM_SUBJECTS).
%
%   Outputs:
%     out - Struct with fields: factor, p_raw, p_GG, eta_p2, n, levels,
%           mu_r, ci_lo_r, ci_hi_r.
%
function out = run_rm_anova(tbl, subjectGroupingVars, factorName, levels, fileID, config)

    % Prepare Fisher z
    r_safe = tbl.Correlation;
    r_safe(r_safe >= 1.0)  = 0.9999999;
    r_safe(r_safe <= -1.0) = -0.9999999;
    tbl.z = atanh(r_safe);

    tbl.Day       = categorical(tbl.Day);
    tbl.Week      = categorical(tbl.Week);
    tbl.TimeOfDay = categorical(tbl.TimeOfDay);

    fprintf('\n[%s] Repeated-Measures ANOVA (subjects = %s)\n', ...
            factorName, strjoin(subjectGroupingVars,' x '));
    fprintf(fileID, '\n[%s] Repeated-Measures ANOVA (subjects = %s)\n', ...
            factorName, strjoin(subjectGroupingVars,' x '));

    % Build mean z per subject x factor level
    S = groupsummary(tbl, [subjectGroupingVars, {factorName}], 'mean', 'z');
    W = unstack(S, 'mean_z', factorName, 'GroupingVariables', subjectGroupingVars);

    % Ensure all level columns exist
    for i = 1:numel(levels)
        if ~ismember(levels{i}, W.Properties.VariableNames)
            W.(levels{i}) = nan(height(W),1);
        end
    end

    % Retain complete subjects
    W = rmmissing(W, 'DataVariables', levels);
    nSubj = height(W);

    if config.LIST_RM_SUBJECTS
        fprintf('  Subjects retained for %s = %d\n', factorName, nSubj);
        fprintf(fileID, '  Subjects retained for %s = %d\n', factorName, nSubj);
        try
            keyStr = evalc('disp(W(:, subjectGroupingVars))');
            fprintf('%s\n', keyStr);
            fprintf(fileID, '%s\n', keyStr);
        catch
        end
    end

    if nSubj < 2
        fprintf('  [SKIP] Not enough subjects.\n');
        fprintf(fileID, '  [SKIP] Not enough subjects.\n');
        out = struct('factor',factorName,'p_raw',NaN,'p_GG',NaN,'eta_p2',NaN, ...
                     'n',nSubj,'levels',{levels}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
        return;
    end

    Y  = W{:, levels};                    % rows: subjects, cols: levels (Fisher z)
    tY = array2table(Y, 'VariableNames', levels);
    withinDesign = table(categorical(levels'), 'VariableNames', {factorName});

    % Assumption checks
    try
        y_stack = Y(:);
        g       = repelem(categorical(levels'), nSubj);
        id      = repelem((1:nSubj)', numel(levels));
        mu_id   = grpstats(y_stack, id, 'mean');
        resid   = y_stack - mu_id(id);

        xs = (resid - mean(resid)) / std(resid);
        [~, pNorm] = adtest(xs);
        pLevene = vartestn(y_stack, g, 'TestType','LeveneAbsolute','Display','off');

        fprintf('    Assumptions: normality p=%.4f (AD); Levene p=%.4f\n', pNorm, pLevene);
        fprintf(fileID, '    Assumptions: normality p=%.4f (AD); Levene p=%.4f\n', pNorm, pLevene);

        if config.SAVE_QQ_PLOTS
            fig = figure('Name', ['QQ Residuals: ' factorName], ...
                         'Visible','off', 'Position',[120,120,700,550]);
            qqplot(resid);
            title(['QQ Plot of Residuals - ' factorName]);
            xlabel('Theoretical Quantiles'); ylabel('Sample Quantiles'); grid on;
            qq_fn = fullfile(config.results_dir, ['qq_residuals_' lower(factorName)]);
            saveas(fig, [qq_fn '.png']);
            saveas(fig, [qq_fn '.svg']);
            savefig(fig, [qq_fn '.fig']);
            close(fig);
        end
    catch MEa
        fprintf('    [WARN] Assumption checks skipped: %s\n', MEa.message);
        fprintf(fileID, '    [WARN] Assumption checks skipped: %s\n', MEa.message);
    end

    try
        % Fit RM model
        formula = sprintf('%s-%s ~ 1', levels{1}, levels{end});
        rm  = fitrm(tY, formula, 'WithinDesign', withinDesign);
        ran = ranova(rm, 'WithinModel', factorName);

        % Locate factor and error rows
        rn       = string(ran.Properties.RowNames);
        termName = "(Intercept):" + string(factorName);
        idx      = find(rn == termName, 1);
        if isempty(idx)
            idx = find(contains(rn, string(factorName)), 1, 'first');
        end
        errIdx = find(contains(rn, "Error(" + string(factorName) + ")"), 1, 'first');
        if isempty(errIdx), errIdx = idx + 1; end
        if isempty(idx)
            error('Could not locate factor row for %s.', factorName);
        end

        % F, p, effect size
        df1   = ran.DF(idx);
        df2   = ran.DF(errIdx);
        Fval  = ran.F(idx);
        p_raw = ran.pValue(idx);

        p_GG = p_raw;
        if ismember('pValueGG', ran.Properties.VariableNames) && ~isnan(ran.pValueGG(idx))
            p_GG = ran.pValueGG(idx);
        end

        SS_eff = ran.SumSq(idx);
        SS_err = ran.SumSq(errIdx);
        eta_p2 = SS_eff / (SS_eff + SS_err);

        % Optional epsilon (log only)
        try
            epsTbl = epsilon(rm);
            eRow = epsTbl(strcmp(string(epsTbl.WithinEffect), factorName), :);
        catch
            eRow = [];
        end

        % RM-ANOVA summary table
        fprintf('\n    RM-ANOVA Summary (%s)\n', factorName);
        fprintf('    %-10s | %4s | %4s | %8s | %10s | %10s | %10s\n', ...
            'Effect','df1','df2','F','p','p(GG)','eta_p^2');
        fprintf('    %-10s | %4d | %4d | %8.3f | %10.4f | %10.4f | %10.3f\n', ...
            factorName, df1, df2, Fval, p_raw, p_GG, eta_p2);

        fprintf(fileID, '\n    RM-ANOVA Summary (%s)\n', factorName);
        fprintf(fileID, ...
            '    %-10s | %4s | %4s | %8s | %10s | %10s | %10s\n', ...
            'Effect','df1','df2','F','p','p(GG)','eta_p^2');
        fprintf(fileID, ...
            '    %-10s | %4d | %4d | %8.3f | %10.4f | %10.4f | %10.3f\n', ...
            factorName, df1, df2, Fval, p_raw, p_GG, eta_p2);

        if ~isempty(eRow)
            fprintf('    Epsilon (GG/HF/LB) = [%.3f / %.3f / %.3f]\n', ...
                eRow.GreenhouseGeisser, eRow.HuynhFeldt, eRow.LowerBound);
            fprintf(fileID, '    Epsilon (GG/HF/LB) = [%.3f / %.3f / %.3f]\n', ...
                eRow.GreenhouseGeisser, eRow.HuynhFeldt, eRow.LowerBound);
        end

        % ----- Level means and 95% CI on z-scale -----
        mu_z    = mean(Y,1);
        sd_z    = std(Y,0,1);
        se_z    = sd_z ./ sqrt(nSubj);
        dfL     = max(nSubj - 1,1);
        tcrit   = tinv(0.975, dfL);
        ci_lo_z = mu_z - tcrit .* se_z;
        ci_hi_z = mu_z + tcrit .* se_z;

        % Back-transform to r-scale
        mu_r    = tanh(mu_z);
        ci_lo_r = tanh(ci_lo_z);
        ci_hi_r = tanh(ci_hi_z);

        fprintf('    Level means (z''''): ');
        for j = 1:numel(levels), fprintf('%s=%.3f ', levels{j}, mu_z(j)); end
        fprintf('\n');

        fprintf('    Level means (r from z): ');
        for j = 1:numel(levels), fprintf('%s=%.4f ', levels{j}, mu_r(j)); end
        fprintf('\n');

        fprintf('    Level 95%% CI (r from z): ');
        for j = 1:numel(levels)
            fprintf('%s=[%.4f,%.4f] ', levels{j}, ci_lo_r(j), ci_hi_r(j));
        end
        fprintf('\n');

        fprintf(fileID, '    Level means (z''''): ');
        for j = 1:numel(levels), fprintf(fileID, '%s=%.3f ', levels{j}, mu_z(j)); end
        fprintf(fileID, '\n');

        fprintf(fileID, '    Level means (r from z): ');
        for j = 1:numel(levels), fprintf(fileID, '%s=%.4f ', levels{j}, mu_r(j)); end
        fprintf(fileID, '\n');

        fprintf(fileID, '    Level 95%% CI (r from z): ');
        for j = 1:numel(levels)
            fprintf(fileID, '%s=[%.4f,%.4f] ', levels{j}, ci_lo_r(j), ci_hi_r(j));
        end
        fprintf(fileID, '\n');

        % Dot-and-whisker plot for this factor
        try
            plot_dotwhisker_r(mu_r, ci_lo_r, ci_hi_r, levels, ...
                config.results_dir, lower(factorName), ...
                ['Temporal reliability: ' factorName]);
        catch MEp
            fprintf('    [WARN] Dot-whisker plot failed for %s: %s\n', factorName, MEp.message);
            fprintf(fileID,'    [WARN] Dot-whisker plot failed for %s: %s\n', factorName, MEp.message);
        end

        out = struct( ...
            'factor',   factorName, ...
            'p_raw',    p_raw, ...
            'p_GG',     p_GG, ...
            'eta_p2',   eta_p2, ...
            'n',        nSubj, ...
            'levels',   {levels}, ...
            'mu_r',     mu_r, ...
            'ci_lo_r',  ci_lo_r, ...
            'ci_hi_r',  ci_hi_r);

    catch MEf
        fprintf('    [ERROR] RM-ANOVA failed for %s: %s\n', factorName, MEf.message);
        fprintf(fileID, '    [ERROR] RM-ANOVA failed for %s: %s\n', factorName, MEf.message);
        out = struct('factor',factorName,'p_raw',NaN,'p_GG',NaN,'eta_p2',NaN, ...
                     'n',nSubj,'levels',{levels}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
    end
end

% run_week_paired - Paired t-test and Wilcoxon signed-rank for WK01 vs WK02.
%
%   Groups rows by (Day, TimeOfDay, Week), computes mean Fisher z per cell,
%   unstacks into paired WK01/WK02 columns, runs a paired t-test and a
%   Wilcoxon signed-rank test on the differences, and saves a dot-whisker
%   plot of week-level means on the r scale.
%
%   Inputs:
%     tbl    - Results table with columns: Correlation, Week, Day, TimeOfDay.
%     fileID - Log file descriptor.
%     config - Config struct (needs results_dir).
%
%   Outputs:
%     out - Struct with fields: factor, p_t, p_wil, df, t_stat, n, levels,
%           mu_r, ci_lo_r, ci_hi_r.
%
function out = run_week_paired(tbl, fileID, config)

    r_safe = tbl.Correlation;
    r_safe(r_safe >= 1.0)  = 0.9999999;
    r_safe(r_safe <= -1.0) = -0.9999999;
    tbl.z = atanh(r_safe);

    tbl.Day       = categorical(tbl.Day);
    tbl.Week      = categorical(tbl.Week);
    tbl.TimeOfDay = categorical(tbl.TimeOfDay);

    fprintf('\n[Week] Paired comparison WK01 vs WK02 (matched Day x TimeOfDay)\n');
    fprintf(fileID, '\n[Week] Paired comparison WK01 vs WK02 (matched Day x TimeOfDay)\n');

    S = groupsummary(tbl, {'Day','TimeOfDay','Week'}, 'mean', 'z');
    W = unstack(S, 'mean_z', 'Week', 'GroupingVariables', {'Day','TimeOfDay'});

    if ~ismember('WK01', W.Properties.VariableNames), W.WK01 = nan(height(W),1); end
    if ~ismember('WK02', W.Properties.VariableNames), W.WK02 = nan(height(W),1); end

    W = rmmissing(W, 'DataVariables', {'WK01','WK02'});
    nPairs = height(W);

    if nPairs < 2
        fprintf('  [SKIP] Not enough matched pairs.\n');
        fprintf(fileID, '  [SKIP] Not enough matched pairs.\n');
        out = struct('factor','Week','p_t',NaN,'p_wil',NaN,'df',NaN,'t_stat',NaN, ...
                     'n',nPairs,'levels',{{'WK01','WK02'}}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
        return;
    end

    x = W.WK01;        % mean z per condition for WK01
    y = W.WK02;        % mean z per condition for WK02
    d = x - y;

    % Normality (AD) on differences
    try
        xs = (d - mean(d)) / std(d);
        [~, pNorm] = adtest(xs);
        fprintf('    Normality (diff) p=%.4f (AD)\n', pNorm);
        fprintf(fileID, '    Normality (diff) p=%.4f (AD)\n', pNorm);
    catch MEa
        fprintf('    [WARN] Normality check skipped: %s\n', MEa.message);
        fprintf(fileID, '    [WARN] Normality check skipped: %s\n', MEa.message);
    end

    try
        [~, p_t, ci_z, stats] = ttest(x, y);
        mean_diff_z = mean(d);
        p_wil = signrank(x, y);

        fprintf('\n    Week Comparison Summary (z-scale)\n');
        fprintf('    %-12s | %4s | %10s | %10s | %25s\n', ...
            'Test','df','t','p','Mean diff z [95% CI]');
        fprintf('    %-12s | %4d | %10.4f | %10.4f | %8.4f [%8.4f, %8.4f]\n', ...
            'paired t', stats.df, stats.tstat, p_t, mean_diff_z, ci_z(1), ci_z(2));
        fprintf('    %-12s | %4s | %10s | %10.4f | %25s\n', ...
            'Wilcoxon','--','--', p_wil, 'n/a');

        fprintf(fileID, '\n    Week Comparison Summary (z-scale)\n');
        fprintf(fileID, '    %-12s | %4s | %10s | %10s | %25s\n', ...
            'Test','df','t','p','Mean diff z [95% CI]');
        fprintf(fileID, '    %-12s | %4d | %10.4f | %10.4f | %8.4f [%8.4f, %8.4f]\n', ...
            'paired t', stats.df, stats.tstat, p_t, mean_diff_z, ci_z(1), ci_z(2));
        fprintf(fileID, '    %-12s | %4s | %10s | %10.4f | %25s\n', ...
            'Wilcoxon','--','--', p_wil, 'n/a');

        % Back-transformed WK means & 95% CI (for plotting)
        levels = {'WK01','WK02'};
        mu_z   = [mean(x), mean(y)];
        sd_z   = [std(x),  std(y)];
        se_z   = sd_z ./ sqrt(nPairs);
        dfL    = max(nPairs - 1, 1);
        tcrit  = tinv(0.975, dfL);
        ci_lo_z = mu_z - tcrit .* se_z;
        ci_hi_z = mu_z + tcrit .* se_z;

        mu_r    = tanh(mu_z);
        ci_lo_r = tanh(ci_lo_z);
        ci_hi_r = tanh(ci_hi_z);

        fprintf('\n    Week means (r from z): WK01=%.4f, WK02=%.4f\n', ...
            mu_r(1), mu_r(2));
        fprintf(fileID, '\n    Week means (r from z): WK01=%.4f, WK02=%.4f\n', ...
            mu_r(1), mu_r(2));

        % Dot-whisker plot for Week comparison
        try
            plot_dotwhisker_r(mu_r, ci_lo_r, ci_hi_r, levels, ...
                config.results_dir, 'week', ...
                'Week-to-week temporal reliability');
        catch MEp2
            fprintf('    [WARN] Dot-whisker plot failed for week comparison: %s\n', MEp2.message);
            fprintf(fileID,'    [WARN] Dot-whisker plot failed for week comparison: %s\n', MEp2.message);
        end

        out = struct( ...
            'factor',   'Week', ...
            'p_t',      p_t, ...
            'p_wil',    p_wil, ...
            'df',       stats.df, ...
            't_stat',   stats.tstat, ...
            'n',        nPairs, ...
            'levels',   {{'WK01','WK02'}}, ...
            'mu_r',     mu_r, ...
            'ci_lo_r',  ci_lo_r, ...
            'ci_hi_r',  ci_hi_r);

    catch MEp
        fprintf('    [ERROR] Paired tests failed: %s\n', MEp.message);
        fprintf(fileID, '    [ERROR] Paired tests failed: %s\n', MEp.message);
        out = struct('factor','Week','p_t',NaN,'p_wil',NaN,'df',NaN,'t_stat',NaN, ...
                     'n',nPairs,'levels',{{'WK01','WK02'}}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
    end
end

% plot_combined_inferential_boxplots - Save three-panel dot-whisker summary figure.
%
%   Lays out three side-by-side dot-and-whisker subplots showing mean Pearson r
%   and 95% CI for (1) time-of-day, (2) day-of-week, and (3) week-to-week factors.
%   Saves the figure as PNG, SVG, and FIG to the results directory.
%
%   Inputs:
%     out_time    - Output struct from run_rm_anova for the TimeOfDay factor.
%     out_day     - Output struct from run_rm_anova for the Day factor.
%     out_week    - Output struct from run_week_paired.
%     results_dir - Path string of the directory where figures are saved.
%
%   Outputs:
%     (none) - Saves inferential_dotwhisker_combined.{png,svg,fig}.
%
function plot_combined_inferential_boxplots(out_time, out_day, out_week, results_dir)

    fig = figure('Name','Inferential_DotWhisker_Combined', ...
                 'Visible','off', ...
                 'Position',[100,100,1350,420]);

    violet_main = [63  81 181] / 255;
    violet_dark = [26  35 126] / 255;

    % --- Subplot 1: Time-of-day ---
    ax1 = subplot(1,3,1);
    if isfield(out_time,'mu_r') && ~isempty(out_time.mu_r)
        draw_dotwhisker_on_ax(ax1, out_time.mu_r, out_time.ci_lo_r, out_time.ci_hi_r, ...
            out_time.levels, 'Time-of-day', true, violet_main, violet_dark);
    else
        axis(ax1,'off');
        text(0.5,0.5,'Time-of-day N/A','HorizontalAlignment','center');
    end

    % --- Subplot 2: Day-of-week ---
    ax2 = subplot(1,3,2);
    if isfield(out_day,'mu_r') && ~isempty(out_day.mu_r)
        draw_dotwhisker_on_ax(ax2, out_day.mu_r, out_day.ci_lo_r, out_day.ci_hi_r, ...
            out_day.levels, 'Day-of-week', false, violet_main, violet_dark);
    else
        axis(ax2,'off');
        text(0.5,0.5,'Day-of-week N/A','HorizontalAlignment','center');
    end

    % --- Subplot 3: Week-to-week ---
    ax3 = subplot(1,3,3);
    if isfield(out_week,'mu_r') && ~isempty(out_week.mu_r)
        draw_dotwhisker_on_ax(ax3, out_week.mu_r, out_week.ci_lo_r, out_week.ci_hi_r, ...
            out_week.levels, 'Week-to-week', false, violet_main, violet_dark);
    else
        axis(ax3,'off');
        text(0.5,0.5,'Week N/A','HorizontalAlignment','center');
    end

    % Save combined figure
    fn = fullfile(results_dir, 'inferential_dotwhisker_combined');
    saveas(fig, [fn '.png']);
    saveas(fig, [fn '.svg']);
    savefig(fig, [fn '.fig']);
    close(fig);

    fprintf('  Saved inferential_dotwhisker_combined.{png,svg,fig}\n');
end

function draw_dotwhisker_on_ax(ax, mu_r, ci_lo_r, ci_hi_r, levels, title_str, showYLabel, ...
                               violet_main, violet_dark)
    axes(ax);
    hold(ax,'on');

    mu_r    = mu_r(:)';
    ci_lo_r = ci_lo_r(:)';
    ci_hi_r = ci_hi_r(:)';
    K       = numel(mu_r);
    x       = 1:K;

    whiskerColor = violet_main;
    dotFaceColor = violet_dark;
    dotEdgeColor = [1 1 1];

    for i = 1:K
        % vertical CI line
        line(ax, [x(i) x(i)], [ci_lo_r(i) ci_hi_r(i)], ...
            'Color', whiskerColor, ...
            'LineWidth', 2.0);
        % mean dot
        plot(ax, x(i), mu_r(i), 'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', dotFaceColor, ...
            'MarkerEdgeColor', dotEdgeColor, ...
            'LineWidth', 1.2);
    end

    if K > 1
        plot(ax, x, mu_r, '-', ...
            'Color', whiskerColor, ...
            'LineWidth', 1.0);
    end

    set(ax, 'XLim', [0.5 K+0.5], ...
            'XTick', x, ...
            'XTickLabel', levels, ...
            'FontName', 'Times', ...
            'FontSize', 9, ...
            'Box', 'off', ...
            'Layer', 'top', ...
            'YGrid', 'on', ...
            'GridColor', [0.90 0.91 0.97], ...
            'GridAlpha', 0.9, ...
            'LineWidth', 0.8);

    if showYLabel
        ylabel(ax, 'Pearson r', 'FontSize', 9);
    else
        ylabel(ax, '');
    end

    % y-limits based on CI, clipped to [0,1]
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

    title(ax, title_str, 'FontWeight','normal', 'FontSize',10);
    hold(ax,'off');
end


function plot_dotwhisker_r(mu_r, ci_lo_r, ci_hi_r, levels, results_dir, tag, title_str)
    % Dot-and-whisker plot for means 95% CIs on r-scale.
    % Uses violet color scheme (consistent with combined plots).

    if nargin < 7 || isempty(title_str)
        title_str = ['Estimates with 95% CI  ' upper(tag)];
    end

    mu_r    = mu_r(:)';
    ci_lo_r = ci_lo_r(:)';
    ci_hi_r = ci_hi_r(:)';
    K       = numel(mu_r);
    x       = 1:K;

    violet_main = [63  81 181] / 255;
    violet_dark = [26  35 126] / 255;
    whiskerColor = violet_main;
    dotFaceColor = violet_dark;
    dotEdgeColor = [1 1 1];

    fig = figure('Name',['DotWhisker_' tag], ...
                 'Visible','off', ...
                 'Position',[120,120,820,560]);
    ax = gca;
    hold(ax,'on');

    for i = 1:K
        line(ax, [x(i) x(i)], [ci_lo_r(i) ci_hi_r(i)], ...
            'Color', whiskerColor, ...
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
            'GridColor', [0.90 0.91 0.97], ...
            'GridAlpha', 0.9, ...
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
    % Use 0.001 step for y-ticks
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function apply_holm_two_tests(outA, outB, fileID)
    fprintf('\n[Holm-Bonferroni across RM tests (GG p-values)]\n');
    fprintf(fileID, '\n[Holm-Bonferroni across RM tests (GG p-values)]\n');

    pvals = [outA.p_GG, outB.p_GG];
    names = {outA.factor, outB.factor};

    [p_sorted, order] = sort(pvals);
    m = numel(pvals);
    any_sig = false;

    for i = 1:m
        idx     = order(i);
        alpha_i = 0.05 / (m - i + 1);
        decision = (p_sorted(i) < alpha_i);
        any_sig  = any_sig | decision;

        fprintf('  %s: p=%.4f vs alpha_Holm=%.4f -> %s\n', ...
            names{idx}, p_sorted(i), alpha_i, ternary(decision,'SIGNIFICANT','ns'));
        fprintf(fileID, '  %s: p=%.4f vs alpha_Holm=%.4f -> %s\n', ...
            names{idx}, p_sorted(i), alpha_i, ternary(decision,'SIGNIFICANT','ns'));

        if ~decision
            fprintf('  (Holm stops here; remaining are non-significant.)\n');
            fprintf(fileID, '  (Holm stops here; remaining are non-significant.)\n');
            break;
        end
    end

    if ~any_sig
        fprintf('  Result: No RM effect survives Holm correction.\n');
        fprintf(fileID, '  Result: No RM effect survives Holm correction.\n');
    end
end

function print_raw_correlations(results_tbl, all_correlations, fileID)
    fprintf('All correlations calculated.\n\n');
    fprintf(fileID, 'All correlations calculated.\n\n');

    fprintf('\nFull Experiment Summary\n');
    fprintf(fileID, '\nFull Experiment Summary\n');

    try
        unique_paths = unique(results_tbl.FilePath);
        fprintf('  Total Unique Traces: %d\n', numel(unique_paths));
        fprintf('  Total Pairwise Comparisons: %d\n', height(results_tbl));
        fprintf(fileID, '  Total Unique Traces: %d\n', numel(unique_paths));
        fprintf(fileID, '  Total Pairwise Comparisons: %d\n', height(results_tbl));
    catch ME
        fprintf('  Summary error: %s\n', ME.message);
        fprintf(fileID, '  Summary error: %s\n', ME.message);
    end

    fprintf('\nAll Correlations (N=%d)\n', numel(all_correlations));
    fprintf(fileID, '\nAll Correlations (N=%d)\n', numel(all_correlations));
    for i = 1:numel(all_correlations)
        fprintf('%0.6f  ', all_correlations(i));
        fprintf(fileID, '%0.6f  ', all_correlations(i));
        if mod(i,5)==0
            fprintf('\n'); fprintf(fileID,'\n');
        end
    end
    fprintf('\n'); fprintf(fileID,'\n');

    fprintf('\nFull Results Table\n');
    fprintf(fileID, '\nFull Results Table\n');
    try
        T_str = evalc('disp(results_tbl)');
        fprintf('%s\n', T_str);
        fprintf(fileID, '%s\n', T_str);
    catch ME2
        fprintf('  Table print error: %s\n', ME2.message);
        fprintf(fileID, '  Table print error: %s\n', ME2.message);
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function plot_box_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val,'z')
        plot_name = 'Z-Transformed Correlation Box Plot';
        y_label   = 'Fisher z';
        suffix    = '_z';
    else
        plot_name = 'Correlation Box Plot';
        y_label   = 'Pearson r';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible','off', ...
                 'Position',[100,100,700,600]);
    ax = gca;

    if ~isempty(data)
        boxplot(ax, data', 'Notch','on','Labels',{'All'});
        margin = (stats.max - stats.min)*0.1;
        if margin <= 0 || isnan(margin), margin = 0.001; end
        ylim(ax, [stats.min - margin, stats.max + margin]);
        ylabel(ax, y_label);
        title(ax, plot_name);
        grid(ax,'on');
    else
        title(ax,'No Data');
        set(ax,'XTick',[],'YTick',[]);
    end

    saveas(fig, fullfile(results_dir, ['boxplot' suffix '.png']));
    savefig(fig, fullfile(results_dir, ['boxplot' suffix '.fig']));
    saveas(fig, fullfile(results_dir, ['boxplot' suffix '.svg']));
    fprintf('Saved boxplot%s.{png,fig}\n', suffix);
    fprintf(fileID,'Saved boxplot%s.{png,fig}\n', suffix);
    close(fig);
end

function plot_hist_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val,'z')
        plot_name = 'Z-Transformed Correlation Histogram';
        x_label   = 'Fisher z';
        suffix    = '_z';
    else
        plot_name = 'Correlation Histogram';
        x_label   = 'Pearson r';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible','off', ...
                 'Position',[100,100,700,500]);
    ax = gca;

    if ~isempty(data)
        histogram(ax, data, 'BinMethod','auto');
        margin = (stats.max - stats.min)*0.1;
        if margin <= 0 || isnan(margin), margin = 0.001; end
        xlim(ax, [stats.min - margin, stats.max + margin]);
    end

    title(ax, [plot_name ' (Zoomed)']);
    xlabel(ax, x_label);
    ylabel(ax, 'Count');
    grid(ax,'on');

    saveas(fig, fullfile(results_dir, ['hist' suffix '.png']));
    savefig(fig, fullfile(results_dir, ['hist' suffix '.fig']));
    saveas(fig, fullfile(results_dir, ['hist' suffix '.svg']));
    fprintf('Saved hist%s.{png,fig}\n', suffix);
    fprintf(fileID,'Saved hist%s.{png,fig}\n', suffix);
    close(fig);
end

function plot_z_score_distributions(z_data, z_stats, results_dir, fileID)
    if isempty(z_data) || ~isfinite(z_stats.std)
        return;
    end

    z_mean = z_stats.mean;
    z_std  = z_stats.std;
    if z_std == 0, z_std = 1e-6; end
    SE   = z_stats.sem;
    df   = max(z_stats.N - 1, 1);
    tcrit= tinv(0.975, df);
    CI   = [z_mean - tcrit*SE, z_mean + tcrit*SE];

    min_val = z_mean - 4*z_std;
    max_val = z_mean + 4*z_std;
    if ~isfinite(min_val), min_val = -3; end
    if ~isfinite(max_val), max_val = 3;  end
    if min_val == max_val, min_val = min_val-1; max_val = max_val+1; end

    x   = linspace(min_val, max_val, 1000);
    pdf = normpdf(x, z_mean, z_std);

    fig = figure('Name','Z-Score Distribution','Visible','off', ...
                 'Position',[100,100,1000,600]);
    ax = gca; hold(ax,'on');
    plot(ax, x, pdf, 'LineWidth',2);
    yl = ylim(ax);
    line(ax, [z_mean z_mean], yl, 'LineStyle','--','LineWidth',1.5);
    line(ax, [CI(1) CI(1)], yl, 'LineWidth',2);
    line(ax, [CI(2) CI(2)], yl, 'LineWidth',2);
    ylim(ax, yl);
    title(ax, 'Normal Approximation of Fisher z');
    xlabel(ax, 'z'); ylabel(ax, 'PDF');
    grid(ax,'on'); hold(ax,'off');

    saveas(fig, fullfile(results_dir, 'z_score_distributions.png'));
    savefig(fig, fullfile(results_dir, 'z_score_distributions.fig'));
    saveas(fig, fullfile(results_dir, 'z_score_distributions.svg'));
    fprintf('Saved z_score_distributions.{png,fig}\n');
    fprintf(fileID,'Saved z_score_distributions.{png,fig}\n');
    close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replay Attack Resilience Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run_replay_attack_analysis - Orchestrate the older-reference stale-replay test.
%
%   Iterates over each configured target ambient trace, calls
%   replay_run_single_target to compare it against all strictly older mains
%   reference windows, pools the resulting correlations, and reports
%   aggregate authentic-vs-replay statistics, a per-target CSV summary, and
%   a representative ENF overlay figure.
%
%   Inputs:
%     config - Configuration struct (needs replay_target_specs,
%              replay_overlay_target_spec, and all ENF extraction fields).
%     fileID - Log file descriptor.
%
%   Outputs:
%     (none) - Saves CSV summary and overlay figure; prints stats to log.
%
function run_replay_attack_analysis(config, fileID)

    replay_results_dir = 'exp_results/replay_attack_resilience';
    if ~exist(replay_results_dir, 'dir'), mkdir(replay_results_dir); end

    target_subfolders = replay_get_target_subfolders(config);
    if isempty(target_subfolders)
        fprintf('  [WARN] Replay: no target folders configured.\n');
        fprintf(fileID, '  [WARN] Replay: no target folders configured.\n');
        return;
    end

    fprintf('\n  Replay: analyzing %d ambient target folders using older reference windows only...\n', numel(target_subfolders));
    fprintf(fileID, '\n  Replay: analyzing %d ambient target folders using older reference windows only...\n', numel(target_subfolders));

    target_results      = repmat(replay_empty_result_struct(), 0, 1);
    pooled_replay_corrs = [];
    authentic_controls  = [];

    for i = 1:numel(target_subfolders)
        target_folder = target_subfolders{i};
        fprintf('\n  Replay target %d/%d: %s\n', i, numel(target_subfolders), target_folder);
        fprintf(fileID, '\n  Replay target %d/%d: %s\n', i, numel(target_subfolders), target_folder);

        res = replay_run_single_target(config, target_folder, fileID);
        target_results(end+1,1) = res; %#ok<AGROW>

        if ~isempty(res.replay_corrs)
            pooled_replay_corrs = [pooled_replay_corrs; res.replay_corrs(:)]; %#ok<AGROW>
        end
        if isfinite(res.authentic_r)
            authentic_controls(end+1,1) = res.authentic_r; %#ok<AGROW>
        end
    end

    replay_stats    = sub_calc_stats(pooled_replay_corrs);
    authentic_stats = sub_calc_stats(authentic_controls);

    fprintf('\n  --- Replay Attack Resilience: Pooled Older-Reference Stale-Replay Statistics ---\n');
    fprintf(fileID, '\n  --- Replay Attack Resilience: Pooled Older-Reference Stale-Replay Statistics ---\n');
    fprintf('  Targets requested               | %10d\n', numel(target_subfolders));
    fprintf('  Targets with valid control r    | %10d\n', numel(authentic_controls));
    fprintf(fileID, '  Targets requested               | %10d\n', numel(target_subfolders));
    fprintf(fileID, '  Targets with valid control r    | %10d\n', numel(authentic_controls));
    print_stats_table(replay_stats, fileID, 'Pooled Older-Reference Replay Correlations (r)');

    fprintf('\n  --- Replay Attack Resilience: Authentic Control Statistics ---\n');
    fprintf(fileID, '\n  --- Replay Attack Resilience: Authentic Control Statistics ---\n');
    print_stats_table(authentic_stats, fileID, 'Authentic Controls (r)');

    try
        summary_tbl = replay_build_target_summary_table(target_results);
        summary_csv = fullfile(replay_results_dir, 'replay_attack_target_summary.csv');
        writetable(summary_tbl, summary_csv);
        fprintf('  Saved %s\n', summary_csv);
        fprintf(fileID, '  Saved %s\n', summary_csv);
    catch MEsum
        fprintf('  [WARN] Replay summary CSV failed: %s\n', MEsum.message);
        fprintf(fileID, '  [WARN] Replay summary CSV failed: %s\n', MEsum.message);
    end

    try
        replay_generate_representative_overlay(config, replay_results_dir, fileID);
    catch MEfig
        fprintf('  [ERROR] Replay overlay plot failed: %s\n', MEfig.message);
        fprintf(fileID, '  [ERROR] Replay overlay plot failed: %s\n', MEfig.message);
    end

    fprintf('\n  STEP 3 (older-reference stale replay) complete.\n');
    fprintf(fileID, '\n  STEP 3 (older-reference stale replay) complete.\n');
end


% replay_run_single_target - Compare one ambient trace against all older reference windows.
%
%   Authenticates the ambient trace against its co-located mains reference
%   (authentic control), then compares it against every mains reference window
%   in the dataset whose time ordinal is strictly less than the target's ordinal.
%   Returns a result struct summarising replay and authentic correlations.
%
%   Inputs:
%     config        - Configuration struct with all ENF extraction parameters
%                     and file naming fields.
%     target_folder - Full path to the leaf folder of the target ambient trace.
%     fileID        - Log file descriptor.
%
%   Outputs:
%     result - Struct with fields: target_folder, leaf_folder, authentic_r,
%              replay_corrs, older_refs_available, n_replay, replay_mean,
%              replay_median, replay_std, replay_min, replay_max.
%
function result = replay_run_single_target(config, target_folder, fileID)
    result = replay_empty_result_struct();
    result.target_folder = string(target_folder);

    [~, leaf_folder] = fileparts(target_folder);
    result.leaf_folder = string(leaf_folder);

    target_meta = replay_parse_folder_metadata(target_folder);
    if ~isfinite(target_meta.ordinal)
        fprintf('  [WARN] Replay: could not parse target folder metadata: %s\n', target_folder);
        fprintf(fileID, '  [WARN] Replay: could not parse target folder metadata: %s\n', target_folder);
        return;
    end

    ambient_wav       = fullfile(target_folder, [config.file_2_name_base '.wav']);
    authentic_ref_wav = fullfile(target_folder, [config.file_1_name_base '.wav']);

    if ~isfile(ambient_wav)
        fprintf('  [WARN] Replay: ambient WAV not found: %s\n', ambient_wav);
        fprintf(fileID, '  [WARN] Replay: ambient WAV not found: %s\n', ambient_wav);
        return;
    end
    if ~isfile(authentic_ref_wav)
        fprintf('  [WARN] Replay: authentic reference WAV not found: %s\n', authentic_ref_wav);
        fprintf(fileID, '  [WARN] Replay: authentic reference WAV not found: %s\n', authentic_ref_wav);
        return;
    end

    try
        result.authentic_r = proc_enf_analysis( ...
            ambient_wav, authentic_ref_wav, ...
            config.nfft, config.frame_size, config.overlap_size, ...
            config.harmonics_arr_1, config.nominal_freq_arr(2), ...
            config.harmonics_arr_2, config.nominal_freq_arr(2), ...
            config.trace_1_freq_est_method, config.trace_1_est_freq, ...
            config.trace_1_freq_est_spec_comb_harmonics, ...
            config.trace_2_freq_est_method, config.trace_2_est_freq, ...
            config.trace_2_freq_est_spec_comb_harmonics, ...
            config.trace_1_plot_title, config.trace_2_plot_title, false);
        fprintf('  Replay: authentic control r=%.6f (%s)\n', result.authentic_r, authentic_ref_wav);
        fprintf(fileID, '  Replay: authentic control r=%.6f (%s)\n', result.authentic_r, authentic_ref_wav);
    catch MEauth
        fprintf('  [WARN] Replay: authentic control failed for %s: %s\n', target_folder, MEauth.message);
        fprintf(fileID, '  [WARN] Replay: authentic control failed for %s: %s\n', target_folder, MEauth.message);
    end

    all_ref_wavs = replay_build_reference_bank(config, leaf_folder);
    older_ref_wavs = cell(0,1);
    for k = 1:numel(all_ref_wavs)
        ref_folder = fileparts(all_ref_wavs{k});
        ref_meta = replay_parse_folder_metadata(ref_folder);
        if isfinite(ref_meta.ordinal) && ref_meta.ordinal < target_meta.ordinal
            older_ref_wavs{end+1,1} = all_ref_wavs{k}; %#ok<AGROW>
        end
    end

    result.older_refs_available = sum(cellfun(@(p) isfile(p), older_ref_wavs));
    fprintf('  Replay: comparing against %d older reference windows only (%d present on disk).\n', ...
        numel(older_ref_wavs), result.older_refs_available);
    fprintf(fileID, '  Replay: comparing against %d older reference windows only (%d present on disk).\n', ...
        numel(older_ref_wavs), result.older_refs_available);

    replay_corrs = [];
    for k = 1:numel(older_ref_wavs)
        ref_wav = older_ref_wavs{k};

        if ~isfile(ref_wav)
            fprintf('  Replay: missing older ref %d/%d - %s\n', k, numel(older_ref_wavs), ref_wav);
            fprintf(fileID, '  Replay: missing older ref %d/%d - %s\n', k, numel(older_ref_wavs), ref_wav);
            continue;
        end

        try
            r_val = proc_enf_analysis( ...
                ambient_wav, ref_wav, ...
                config.nfft, config.frame_size, config.overlap_size, ...
                config.harmonics_arr_1, config.nominal_freq_arr(2), ...
                config.harmonics_arr_2, config.nominal_freq_arr(2), ...
                config.trace_1_freq_est_method, config.trace_1_est_freq, ...
                config.trace_1_freq_est_spec_comb_harmonics, ...
                config.trace_2_freq_est_method, config.trace_2_est_freq, ...
                config.trace_2_freq_est_spec_comb_harmonics, ...
                config.trace_1_plot_title, config.trace_2_plot_title, false);
            replay_corrs(end+1,1) = r_val; %#ok<AGROW>
            fprintf('  Replay: older ref %d/%d  r=%.6f  %s\n', k, numel(older_ref_wavs), r_val, ref_wav);
            fprintf(fileID, '  Replay: older ref %d/%d  r=%.6f  %s\n', k, numel(older_ref_wavs), r_val, ref_wav);
        catch MEk
            fprintf('  Replay: error older ref %d/%d - %s\n', k, numel(older_ref_wavs), MEk.message);
            fprintf(fileID, '  Replay: error older ref %d/%d - %s\n', k, numel(older_ref_wavs), MEk.message);
        end
    end

    replay_corrs = replay_corrs(isfinite(replay_corrs));
    replay_stats = sub_calc_stats(replay_corrs);

    result.replay_corrs  = replay_corrs;
    result.n_replay      = replay_stats.N;
    result.replay_mean   = replay_stats.mean;
    result.replay_median = replay_stats.median;
    result.replay_std    = replay_stats.std;
    result.replay_min    = replay_stats.min;
    result.replay_max    = replay_stats.max;

    if result.n_replay > 0
        fprintf('  Replay target summary: older refs=%d, N=%d, mean=%.6f, median=%.6f, std=%.6f\n', ...
            result.older_refs_available, result.n_replay, result.replay_mean, result.replay_median, result.replay_std);
        fprintf(fileID, '  Replay target summary: older refs=%d, N=%d, mean=%.6f, median=%.6f, std=%.6f\n', ...
            result.older_refs_available, result.n_replay, result.replay_mean, result.replay_median, result.replay_std);
    elseif result.older_refs_available == 0
        fprintf('  [WARN] Replay: no older reference windows available for %s\n', target_folder);
        fprintf(fileID, '  [WARN] Replay: no older reference windows available for %s\n', target_folder);
    else
        fprintf('  [WARN] Replay: no valid older-reference stale-replay correlations for %s\n', target_folder);
        fprintf(fileID, '  [WARN] Replay: no valid older-reference stale-replay correlations for %s\n', target_folder);
    end
end


function target_subfolders = replay_get_target_subfolders(config)
    target_subfolders = cell(0,1);
    if ~isfield(config, 'replay_target_specs') || isempty(config.replay_target_specs)
        return;
    end

    target_specs = config.replay_target_specs;
    for i = 1:size(target_specs, 1)
        target_subfolders{end+1,1} = replay_spec_to_subfolder(config, target_specs(i,:)); %#ok<AGROW>
    end
end


function subfolder = replay_spec_to_subfolder(config, spec)
    week      = char(spec{1});
    day_type  = char(spec{2});
    day_token = upper(char(spec{3}));
    time_slot = char(spec{4});
    leaf      = char(spec{5});

    if strcmp(day_token, 'THU') || strcmp(day_token, 'THUR')
        day_token = replay_get_thursday_token(config);
    end

    subfolder = fullfile(config.baseDir, week, day_type, day_token, time_slot, leaf);
end


function all_ref_wavs = replay_build_reference_bank(config, leaf_folder)
    leaf_folder = char(leaf_folder);
    ref_name   = [config.file_1_name_base '.wav'];
    weeks      = {'WK01', 'WK02'};
    thu_token  = replay_get_thursday_token(config);
    day_paths  = {'WEND/SAT', 'WEND/SUN', 'WDAY/WED', ['WDAY/' thu_token]};
    time_slots = {'EMRN', 'MORN', 'AFTN', 'EVEN'};

    all_ref_wavs = cell(0,1);
    for wi = 1:numel(weeks)
        for di = 1:numel(day_paths)
            for ti = 1:numel(time_slots)
                p = fullfile(config.baseDir, weeks{wi}, day_paths{di}, time_slots{ti}, leaf_folder, ref_name);
                all_ref_wavs{end+1,1} = p; %#ok<AGROW>
            end
        end
    end
end


function ordinal = replay_get_time_ordinal(week_str, ~, day_str, slot_str)
    week_str = upper(char(week_str));
    day_str  = upper(char(day_str));
    slot_str = upper(char(slot_str));

    week_idx = NaN;
    switch week_str
        case 'WK01'
            week_idx = 1;
        case 'WK02'
            week_idx = 2;
    end

    day_idx = NaN;
    switch day_str
        case 'WED'
            day_idx = 1;
        case {'THU', 'THUR'}
            day_idx = 2;
        case 'SAT'
            day_idx = 3;
        case 'SUN'
            day_idx = 4;
    end

    slot_idx = NaN;
    switch slot_str
        case 'EMRN'
            slot_idx = 1;
        case 'MORN'
            slot_idx = 2;
        case 'AFTN'
            slot_idx = 3;
        case 'EVEN'
            slot_idx = 4;
    end

    if ~isfinite(week_idx) || ~isfinite(day_idx) || ~isfinite(slot_idx)
        ordinal = NaN;
        return;
    end

    ordinal = ((week_idx - 1) * 4 + (day_idx - 1)) * 4 + slot_idx;
end


function meta = replay_parse_folder_metadata(folder_path)
    tokens = upper(regexp(char(folder_path), '[^\\/]+', 'match'));

    week = 'NA';
    if any(strcmp(tokens, 'WK01'))
        week = 'WK01';
    elseif any(strcmp(tokens, 'WK02'))
        week = 'WK02';
    end

    day_type = 'NA';
    if any(strcmp(tokens, 'WDAY'))
        day_type = 'WDAY';
    elseif any(strcmp(tokens, 'WEND'))
        day_type = 'WEND';
    end

    day = 'NA';
    if any(strcmp(tokens, 'WED'))
        day = 'WED';
    elseif any(strcmp(tokens, 'THU')) || any(strcmp(tokens, 'THUR'))
        day = 'THU';
    elseif any(strcmp(tokens, 'SAT'))
        day = 'SAT';
    elseif any(strcmp(tokens, 'SUN'))
        day = 'SUN';
    end

    slot = 'NA';
    if any(strcmp(tokens, 'EMRN'))
        slot = 'EMRN';
    elseif any(strcmp(tokens, 'MORN'))
        slot = 'MORN';
    elseif any(strcmp(tokens, 'AFTN'))
        slot = 'AFTN';
    elseif any(strcmp(tokens, 'EVEN'))
        slot = 'EVEN';
    end

    meta = struct( ...
        'week', string(week), ...
        'day_type', string(day_type), ...
        'day', string(day), ...
        'slot', string(slot), ...
        'ordinal', replay_get_time_ordinal(week, day_type, day, slot));
end


function token = replay_get_thursday_token(config)
    has_thu = exist(fullfile(config.baseDir, 'WK01', 'WDAY', 'THU'), 'dir') || ...
              exist(fullfile(config.baseDir, 'WK02', 'WDAY', 'THU'), 'dir');
    has_thur = exist(fullfile(config.baseDir, 'WK01', 'WDAY', 'THUR'), 'dir') || ...
               exist(fullfile(config.baseDir, 'WK02', 'WDAY', 'THUR'), 'dir');

    token = 'THU';
    if ~has_thu && has_thur
        token = 'THUR';
    end
end


function result = replay_empty_result_struct()
    result = struct( ...
        'target_folder', "", ...
        'leaf_folder', "", ...
        'authentic_r', NaN, ...
        'replay_corrs', [], ...
        'older_refs_available', 0, ...
        'n_replay', 0, ...
        'replay_mean', NaN, ...
        'replay_median', NaN, ...
        'replay_std', NaN, ...
        'replay_min', NaN, ...
        'replay_max', NaN);
end

function summary_tbl = replay_build_target_summary_table(target_results)
    n = numel(target_results);
    summary_tbl = table('Size', [n, 10], ...
                        'VariableTypes', {'string','string','double','double','double','double','double','double','double','double'}, ...
                        'VariableNames', {'TargetFolder','LeafFolder','AuthenticR','OlderRefsAvailable','NReplay','ReplayMean','ReplayMedian','ReplayStd','ReplayMin','ReplayMax'});

    for i = 1:n
        summary_tbl.TargetFolder(i)        = string(target_results(i).target_folder);
        summary_tbl.LeafFolder(i)          = string(target_results(i).leaf_folder);
        summary_tbl.AuthenticR(i)          = target_results(i).authentic_r;
        summary_tbl.OlderRefsAvailable(i) = target_results(i).older_refs_available;
        summary_tbl.NReplay(i)             = target_results(i).n_replay;
        summary_tbl.ReplayMean(i)          = target_results(i).replay_mean;
        summary_tbl.ReplayMedian(i)        = target_results(i).replay_median;
        summary_tbl.ReplayStd(i)           = target_results(i).replay_std;
        summary_tbl.ReplayMin(i)           = target_results(i).replay_min;
        summary_tbl.ReplayMax(i)           = target_results(i).replay_max;
    end
end

function replay_generate_representative_overlay(config, replay_results_dir, fileID)
    if isfield(config, 'replay_overlay_target_spec') && ~isempty(config.replay_overlay_target_spec)
        ambient_subfolder = replay_spec_to_subfolder(config, config.replay_overlay_target_spec);
    else
        target_subfolders = replay_get_target_subfolders(config);
        if isempty(target_subfolders)
            fprintf('  [WARN] Replay overlay skipped: no target folders configured.\n');
            fprintf(fileID, '  [WARN] Replay overlay skipped: no target folders configured.\n');
            return;
        end
        ambient_subfolder = target_subfolders{1};
    end

    target_meta = replay_parse_folder_metadata(ambient_subfolder);
    if ~isfinite(target_meta.ordinal)
        fprintf('  [WARN] Replay overlay skipped: could not parse target time metadata for %s\n', ambient_subfolder);
        fprintf(fileID, '  [WARN] Replay overlay skipped: could not parse target time metadata for %s\n', ambient_subfolder);
        return;
    end

    ambient_wav       = fullfile(ambient_subfolder, [config.file_2_name_base '.wav']);
    authentic_ref_wav = fullfile(ambient_subfolder, [config.file_1_name_base '.wav']);
    [~, leaf_folder]  = fileparts(ambient_subfolder);

    if ~isfile(ambient_wav)
        fprintf('  [WARN] Replay overlay: ambient WAV not found: %s\n', ambient_wav);
        fprintf(fileID, '  [WARN] Replay overlay: ambient WAV not found: %s\n', ambient_wav);
        return;
    end
    if ~isfile(authentic_ref_wav)
        fprintf('  [WARN] Replay overlay: authentic reference WAV not found: %s\n', authentic_ref_wav);
        fprintf(fileID, '  [WARN] Replay overlay: authentic reference WAV not found: %s\n', authentic_ref_wav);
        return;
    end

    try
        [r_auth, ~, ~, ~, enf_ambient, enf_authentic] = proc_enf_analysis( ...
            ambient_wav, authentic_ref_wav, ...
            config.nfft, config.frame_size, config.overlap_size, ...
            config.harmonics_arr_1, config.nominal_freq_arr(2), ...
            config.harmonics_arr_2, config.nominal_freq_arr(2), ...
            config.trace_1_freq_est_method, config.trace_1_est_freq, ...
            config.trace_1_freq_est_spec_comb_harmonics, ...
            config.trace_2_freq_est_method, config.trace_2_est_freq, ...
            config.trace_2_freq_est_spec_comb_harmonics, ...
            config.trace_1_plot_title, config.trace_2_plot_title, false);
    catch MEauth
        fprintf('  [ERROR] Replay overlay: authentic ENF extraction failed: %s\n', MEauth.message);
        fprintf(fileID, '  [ERROR] Replay overlay: authentic ENF extraction failed: %s\n', MEauth.message);
        return;
    end

    sel_subfolders = replay_get_overlay_reference_subfolders(config, leaf_folder);
    sel_short_labels = {
        'WK02 SUN EMRN (same day, earlier)';
        'WK02 THU AFTN';
        'WK01 SUN AFTN (prev week, same day/time)';
        'WK01 WED AFTN';
    };

    sel_enf    = {};
    sel_corrs  = [];
    sel_labels = {};

    for m = 1:numel(sel_subfolders)
        ref_meta = replay_parse_folder_metadata(sel_subfolders{m});
        if ~isfinite(ref_meta.ordinal)
            fprintf('  [WARN] Replay plot: could not parse time metadata for %s\n', sel_subfolders{m});
            fprintf(fileID, '  [WARN] Replay plot: could not parse time metadata for %s\n', sel_subfolders{m});
            continue;
        end
        if ref_meta.ordinal >= target_meta.ordinal
            fprintf('  [WARN] Replay plot: skipping reference that is not older than the target %s for label %s\n', ...
                sel_subfolders{m}, sel_short_labels{m});
            fprintf(fileID, '  [WARN] Replay plot: skipping reference that is not older than the target %s for label %s\n', ...
                sel_subfolders{m}, sel_short_labels{m});
            continue;
        end
        sel_ref_wav = fullfile(sel_subfolders{m}, [config.file_1_name_base '.wav']);
        if ~isfile(sel_ref_wav)
            fprintf('  [WARN] Replay plot: missing %s\n', sel_ref_wav);
            fprintf(fileID, '  [WARN] Replay plot: missing %s\n', sel_ref_wav);
            continue;
        end
        try
            [r_m, ~, ~, ~, ~, enf_m] = proc_enf_analysis( ...
                ambient_wav, sel_ref_wav, ...
                config.nfft, config.frame_size, config.overlap_size, ...
                config.harmonics_arr_1, config.nominal_freq_arr(2), ...
                config.harmonics_arr_2, config.nominal_freq_arr(2), ...
                config.trace_1_freq_est_method, config.trace_1_est_freq, ...
                config.trace_1_freq_est_spec_comb_harmonics, ...
                config.trace_2_freq_est_method, config.trace_2_est_freq, ...
                config.trace_2_freq_est_spec_comb_harmonics, ...
                config.trace_1_plot_title, config.trace_2_plot_title, false);
            sel_enf{end+1}    = enf_m; %#ok<AGROW>
            sel_corrs(end+1)  = r_m; %#ok<AGROW>
            sel_labels{end+1} = sel_short_labels{m}; %#ok<AGROW>
            fprintf('  Replay plot: %s  r=%.6f\n', sel_short_labels{m}, r_m);
            fprintf(fileID, '  Replay plot: %s  r=%.6f\n', sel_short_labels{m}, r_m);
        catch MEsel
            fprintf('  [WARN] Replay plot ENF failed for %s: %s\n', ...
                sel_short_labels{m}, MEsel.message);
            fprintf(fileID, '  [WARN] Replay plot ENF failed for %s: %s\n', ...
                sel_short_labels{m}, MEsel.message);
        end
    end

    traces.enf_ambient      = enf_ambient;
    traces.enf_authentic    = enf_authentic;
    traces.r_authentic      = r_auth;
    traces.sel_enf          = sel_enf;
    traces.sel_corrs        = sel_corrs;
    traces.sel_labels       = sel_labels;
    traces.frame_size_samp  = config.frame_size;
    traces.overlap_samp     = config.overlap_size;

    plot_replay_attack_enf_overlay(traces, replay_results_dir, fileID);
end


function sel_subfolders = replay_get_overlay_reference_subfolders(config, leaf_folder)
    thu_token = replay_get_thursday_token(config);
    sel_subfolders = {
        fullfile(config.baseDir, 'WK02', 'WEND', 'SUN', 'EMRN', leaf_folder), ...
        fullfile(config.baseDir, 'WK02', 'WDAY', thu_token, 'AFTN', leaf_folder), ...
        fullfile(config.baseDir, 'WK01', 'WEND', 'SUN', 'AFTN', leaf_folder), ...
        fullfile(config.baseDir, 'WK01', 'WDAY', 'WED', 'AFTN', leaf_folder), ...
    };
end

function plot_replay_attack_enf_overlay(traces, results_dir, fileID)

    if nargin < 3, fileID = -1; end

    % --- Hop time (seconds per ENF frame) ---
    % frame_size_samp and overlap_samp are in samples at 1000 Hz nominal rate.
    hop_sec = (traces.frame_size_samp - traces.overlap_samp) / 1000;

    % --- Trim all traces to the shortest common length ---
    N = numel(traces.enf_ambient);
    if ~isempty(traces.enf_authentic)
        N = min(N, numel(traces.enf_authentic));
    end
    for m = 1:numel(traces.sel_enf)
        if ~isempty(traces.sel_enf{m})
            N = min(N, numel(traces.sel_enf{m}));
        end
    end
    if N < 2
        fprintf('  [WARN] Replay plot: insufficient ENF frames (%d). Skipping.\n', N);
        if fileID > 0
            fprintf(fileID,'  [WARN] Replay plot: insufficient ENF frames (%d). Skipping.\n', N);
        end
        return;
    end

    t        = (0:N-1) * hop_sec;          % time axis in seconds
    enf_amb  = traces.enf_ambient(1:N);
    enf_auth = traces.enf_authentic(1:N);

    % --- Replay line specifications (color, style, width) ---
    replay_specs = {
        [0.80 0.00 0.80], '-.', 1.6;   % magenta  dash-dot
        [0.00 0.50 0.00], ':',  2.0;   % dark green dotted
        [0.50 0.50 0.50], '--', 1.5;   % gray     dashed
        [0.90 0.48 0.00], '-.',  1.6;  % orange   dash-dot
    };

    % --- Create figure ---
    fig = figure('Name', 'Replay_Attack_ENF_Overlay', ...
                 'Visible', 'off', ...
                 'Units',   'pixels', ...
                 'Position', [100, 100, 1480, 720]);

    % Leave room on the right for the external legend
    ax = axes('Parent', fig, ...
              'Units',  'normalized', ...
              'Position', [0.07, 0.12, 0.60, 0.80]);
    hold(ax, 'on');

    % Plot replay traces first (bottom layer, so authentic and target are on top)
    n_replay = min(numel(traces.sel_enf), size(replay_specs, 1));
    leg_h    = gobjects(2 + n_replay, 1);
    leg_lbl  = cell(2 + n_replay, 1);

    for m = 1:n_replay
        if isempty(traces.sel_enf{m}), continue; end
        enf_m = traces.sel_enf{m}(1:N);
        col   = replay_specs{m, 1};
        ls    = replay_specs{m, 2};
        lw    = replay_specs{m, 3};
        leg_h(2 + m) = plot(ax, t, enf_m, ls, ...
            'Color', col, 'LineWidth', lw);
        leg_lbl{2 + m} = sprintf('Replay: %s ($r = %.3f$)', ...
            traces.sel_labels{m}, traces.sel_corrs(m));
    end

    % Authentic reference (second layer)
    leg_h(2) = plot(ax, t, enf_auth, '--', ...
        'Color', [0.12 0.47 0.71], 'LineWidth', 1.8);
    leg_lbl{2} = sprintf('Authentic Match ($r = %.3f$)', traces.r_authentic);

    % Ambient target on top
    leg_h(1) = plot(ax, t, enf_amb, '-', ...
        'Color', [0.85 0.00 0.00], 'LineWidth', 2.5);
    leg_lbl{1} = 'Ambient Target (Device Trace)';

    hold(ax, 'off');

    % --- Axis formatting ---
    set(ax, ...
        'FontName',  'Times New Roman', ...
        'FontSize',  10, ...
        'Box',       'off', ...
        'LineWidth', 0.8, ...
        'Layer',     'top', ...
        'XGrid',     'on', ...
        'YGrid',     'on', ...
        'GridColor', [0.78 0.78 0.78], ...
        'GridAlpha', 0.75, ...
        'TickDir',   'out', ...
        'XMinorGrid','off', ...
        'YMinorGrid','off');

    xlabel(ax, 'Time (s)', ...
        'FontName', 'Times New Roman', 'FontSize', 11);
    ylabel(ax, 'Instantaneous ENF Value (Hz)', ...
        'FontName', 'Times New Roman', 'FontSize', 11);
    title(ax, 'ENF Replay Attack: Intra-Correlation Resilience Test', ...
        'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'normal');

    xlim(ax, [0, t(end)]);

    % Choose y-limits from the plotted ENF values so peaks are not clipped.
    yvals = [enf_amb(:); enf_auth(:)];
    for m = 1:n_replay
        if ~isempty(traces.sel_enf{m})
            yvals = [yvals; traces.sel_enf{m}(1:N)]; %#ok<AGROW>
        end
    end
    yvals = yvals(isfinite(yvals));
    if isempty(yvals)
        ymin = 59.97;
        ymax = 60.01;
        step = 0.01;
    else
        yrange = max(yvals) - min(yvals);
        pad    = max(0.001, 0.08 * yrange);
        if yrange <= 0.02
            step = 0.002;
        elseif yrange <= 0.05
            step = 0.005;
        else
            step = 0.01;
        end
        ymin = floor((min(yvals) - pad) / step) * step;
        ymax = ceil((max(yvals) + pad) / step) * step;
        if ymin >= ymax
            ymin = ymin - step;
            ymax = ymax + step;
        end
    end
    ylim(ax, [ymin ymax]);
    set(ax, 'YTick', ymin:step:ymax);

    % --- Legend: external, right of axes ---
    valid_mask = arrayfun(@(h) isvalid(h) && ...
        ~isa(h, 'matlab.graphics.GraphicsPlaceholder'), leg_h);
    leg = legend(ax, leg_h(valid_mask), leg_lbl(valid_mask), ...
        'Location',    'eastoutside', ...
        'Interpreter', 'latex', ...
        'FontName',    'Times New Roman', ...
        'FontSize',    9, ...
        'Box',         'on', ...
        'EdgeColor',   [0.75 0.75 0.75]);
    leg.Position(1) = 0.71;   % nudge legend to the right margin

    % --- Save as SVG only (no PNG / FIG) ---
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    fn = fullfile(results_dir, 'replay_attack_enf_overlay');
    saveas(fig, [fn '.svg']);
    close(fig);

    fprintf('  Saved %s.svg\n', fn);
    if fileID > 0
        fprintf(fileID, '  Saved %s.svg\n', fn);
    end
end



