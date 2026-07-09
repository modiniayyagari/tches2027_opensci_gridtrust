%% ENF Pair-Wise Anti-Forensics Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Evaluates three anti-forensics attack variants (A0_2, A1_2, A2_2)
%   applied to one selected authentic AC mains power trace from the artifact
%   dataset. The script generates forged trace variants and then runs two
%   detection techniques to assess how well the forgeries can be identified.
%
% Analysis Steps:
%   1. Load and resample the selected authentic trace from the artifact input tree
%      to the configured analysis sampling rate (ANALYSIS_FS_HZ).
%   2. Apply enabled attack techniques in cascade to generate forged traces:
%        - A0_2: ENF narrow-band bandstop removal (cascade FIR, 2 passes).
%        - A1_2: A0_2 base plus narrow-band noise fill matched to the
%                background spectrum of the authentic trace.
%        - A2_2: A1_2 base plus peak-magnitude-matched FM synthetic ENF
%                component injected at the fundamental frequency.
%   3. Print a consolidated FFT harmonic attenuation comparison report for
%      the authentic vs. each enabled forged trace variant.
%   4. Export per-attack FFT overlay and spectrogram comparison plots.
%   5. Run enabled detection techniques on the forged traces:
%        - D0: Extract ENF signatures from authentic and forged traces and
%              report Pearson correlation with overlay plot.
%        - D1: Compare ENF extracted at the fundamental against ENF extracted
%              at higher harmonics (inter-harmonic consistency check).
%
% Inputs:
%   exp_inputs/MULTI/US_60/OCT/MON/T02/ - Default leaf folder containing
%     the AC mains reference trace mains_pow_trace_ac_egrid_citya_lab.wav.
%     Update selected_leaf_dir and file_1_name in Section 1 to use another trace.
%
% Outputs:
%   exp_results/af_analysis_pairwise/ - Per-attack FFT overlay plots, spectrogram
%     comparisons, ENF overlay plots (D0), and harmonic consistency plots (D1).
%   exp_logs/<script_name>_log.txt    - Full diary log of all console output.
%
% Dependencies:
%   proc_enf_analysis - Pre-compiled ENF extraction and correlation function.
%                       Must be on the MATLAB path before running this script.
%
% Usage:
%   Run from the exp_scripts/ directory. The artifact data root is resolved
%   automatically relative to this script file's location. Update the three
%   variables in Section 1 to point to a different input trace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1) Configuration: Input Trace and STFT Parameters
% Specify the input leaf folder and file name, set the nominal grid frequency,
% and configure the STFT window and overlap parameters used for ENF extraction.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input trace file information.
% This script analyzes one selected authentic trace from one leaf folder.
script_dir = fileparts(mfilename('fullpath'));
artifact_root = fileparts(script_dir);
selected_leaf_dir = fullfile('exp_inputs', 'MULTI', 'US_60', 'OCT', 'MON', 'T02');
file_path = strrep(fullfile(artifact_root, selected_leaf_dir), '\', '/');
file_1_name = "mains_pow_trace_ac_egrid_citya_lab";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram Computation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental power grid frequency settings
nominal_freq_arr      = [50 60];

% Match this value to the selected leaf folder:
%   US_60 -> nominal_freq_arr(2) = 60 Hz
%   DE_50 -> nominal_freq_arr(1) = 50 Hz
nominal_freq_1        = nominal_freq_arr(2);        
harmonics_arr_1       = (1:7)*nominal_freq_1;

% STFT compute param settings
frame_size_arr      = (1:12)*1000;
frame_size          = frame_size_arr(8);                %8000ms window
overlap_size_arr    = 0:0.1:0.9;
overlap_size        = overlap_size_arr(1)*frame_size;   %non-overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Plotting Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the same colormap across generated figures.
set(0,'DefaultFigureColormap', jet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2) Configuration: Anti-Forensics Attack and Detection Settings
% Configure the AF_CFG struct that controls all A0/A1/A2 attack parameters
% and D0/D1 detection settings, including filter bandwidth, FIR order,
% noise-fill statistics, synthetic ENF modulation, and plot export paths.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anti-Forensics (Section 3 A0/A1/A2) Demo Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This block applies anti-forensics operations to the authentic trace (mapped to
% `file_1_name`).
%
% Attack techniques (cascade variants only):
%   A0_2 ENF removal narrow cascade
%   A1_2 ENF noise fill narrow cascade
%   A2_2 ENF embedding narrow cascade

AF_CFG = struct();

% -- Shared -------------------------------------------------------------------
AF_CFG.ANALYSIS_FS_HZ  = 1050;
AF_CFG.NOMINAL_FREQ_HZ = nominal_freq_1;

% Enabled attacks and detectors (single source of truth)
AF_CFG.ATTACK_TECHNIQUES = { ...
    'A0_2 ENF removal narrow cascade', ...
    'A1_2 ENF noise fill narrow cascade', ...
    'A2_2 ENF embedding narrow cascade', ...
};
AF_CFG.DETECTION_TECHNIQUES = { ...
    'D0 enf correlation check', ...
    'D1 inter harmonic consistency check', ...
};

% Harmonic multipliers applied by A0/A1/A2 (fundamental only = [1])
AF_CFG.ATTACK_HARMONIC_MULTS = 1;

% -- A0: ENF Bandstop Removal -------------------------------------------------
AF_CFG.A0_ENF_HALF_BW_HZ          = 1;
AF_CFG.A0_TRANSITION_BW_HZ        = 8;
AF_CFG.A0_FIR_ORDER               = 350;       % ~paper-equivalent at fs=1050 Hz
AF_CFG.A0_FIR_WEIGHTS             = [1 1 1];   % [left stop, pass, right stop]
AF_CFG.A0_CASCADE_PASSES          = 2;
AF_CFG.A0_REPORT_FILTER_RESP_NFFT = 262144;    % dense grid for attenuation report

% -- A1: ENF Noise Fill -------------------------------------------------------
AF_CFG.A1_NOISE_NEIGHBOR_BW_HZ       = 2 * AF_CFG.A0_ENF_HALF_BW_HZ;
AF_CFG.A1_NOISE_MATCH_STAT           = 'median'; % {'median','mean'}
AF_CFG.A1_RANDOM_SEED                = 42;
AF_CFG.A1_USE_EXACT_A0_VARIANT_BASE  = true;

% -- A2: Synthetic ENF Embedding ----------------------------------------------
AF_CFG.A2_SYN_ENF_DEV_HZ    = 0.15;
AF_CFG.A2_SYN_MOD_FREQS_HZ  = [1/120, 1/45];
AF_CFG.A2_SYN_MOD_AMPS_HZ   = [0.08, 0.04];
AF_CFG.A2_RANDOM_SEED        = 42;
AF_CFG.A2_PEAK_MATCH_FREQ_HZ = AF_CFG.NOMINAL_FREQ_HZ;
AF_CFG.A2_PEAK_MATCH_BW_HZ   = AF_CFG.A0_ENF_HALF_BW_HZ;

% -- D0: ENF Correlation Check ------------------------------------------------
AF_CFG.D0_ENF_FREQ_HZ        = nominal_freq_1;
AF_CFG.D0_FRAME_SIZE_SAMPLES = frame_size;
AF_CFG.D0_OVERLAP_SAMPLES    = overlap_size;

% -- D1: Inter-Harmonic Consistency Check -------------------------------------
AF_CFG.D1_HARMONIC_MULTS = [1, 3, 5, 7];
AF_CFG.D1_PLOT_A0_2      = false;
AF_CFG.D1_PLOT_A1        = false;
AF_CFG.D1_PLOT_A2        = true;

% -- Logging / Plot Export ----------------------------------------------------
AF_CFG.LOG_DIR             = fullfile('exp_logs');
AF_CFG.PLOT_EXPORT_DIR     = fullfile('exp_results', 'af_analysis_pairwise');
AF_CFG.PLOT_DISPLAY_ENABLE = false;

% -- FFT Report Settings ------------------------------------------------------
AF_CFG.FFT_HARMONIC_TEXT_FREQS_HZ       = harmonics_arr_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3) Setup: File Paths, Logging, and Plot Export Directory
% Build the full WAV path from the configured leaf folder and file name, then
% open the diary log file and create the plot export directory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input file paths, log file and saved .png file config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_file_1_path_wav = string(fullfile(file_path, file_1_name + ".wav"));

af_script_name = string(mfilename);
if strlength(af_script_name) == 0
    af_script_name = "enf_analysis_top_pair_wise_af";
end
[af_script_dir, ~, ~] = fileparts(mfilename('fullpath'));
if isempty(af_script_dir)
    af_script_dir = pwd;
end
[af_log_file_path, af_plot_export_dir, af_prev_fig_visibility] = ...
    af_prepare_run_logging_and_plot_export(AF_CFG, af_script_name, af_script_dir);
AF_CFG.RUNTIME_PLOT_EXPORT_DIR = af_plot_export_dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 4) Main Analysis: Run Anti-Forensics Demo
% Call af_run_anti_forensics_a0_a1_a2_demo to execute all enabled attack and
% detection techniques on the configured authentic trace. Finalizes the log
% and restores figure visibility on completion or on error.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENF Anti-forensics Analysis (Start)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    fprintf('\nINFO: Running ENF Anti-forensics Analysis\n');
    fprintf('INFO: Authentic trace = %s\n', full_file_1_path_wav);

    af_run_anti_forensics_a0_a1_a2_demo(full_file_1_path_wav, AF_CFG);
    fprintf('\nINFO: ENF Anti-forensics Analysis Complete.\n');
    fprintf('INFO: Log file saved to %s\n', af_log_file_path);
catch ME
    fprintf(2, '\nERROR: %s\n', ME.message);
    fprintf(2, 'ERROR: Stack trace follows in this log session.\n');
    for iE = 1:numel(ME.stack)
        fprintf(2, 'ERROR:   at %s (line %d)\n', ME.stack(iE).name, ME.stack(iE).line);
    end
    af_finalize_run_logging_and_plot_export(af_prev_fig_visibility);
    rethrow(ME);
end
af_finalize_run_logging_and_plot_export(af_prev_fig_visibility);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5) Local Functions: Core Analysis and Helpers
% All local functions used by this script are defined below. The first function
% af_run_anti_forensics_a0_a1_a2_demo is the main analysis entry point; all
% subsequent functions are internal helpers for attacks, detectors, and plotting.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENF Anti-forensics Analysis (Main)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% af_run_anti_forensics_a0_a1_a2_demo - Execute all enabled A0/A1/A2 attack and
%   D0/D1 detection techniques on the authentic trace.
%
%   Loads the authentic trace, runs the enabled cascade attack variants to produce
%   forged traces, prints the FFT harmonic attenuation report, exports per-attack
%   plots, then runs the enabled detectors (D0 ENF correlation, D1 inter-harmonic
%   consistency) on all forged variants.
%
%   Inputs:
%     authentic_wav_path - Full path to the authentic AC mains WAV file.
%     cfg                - AF_CFG struct with all attack and detection parameters.
%
%   Outputs:
%     (none) - All results printed to console (diary log) and saved as plot files.
%
function af_run_anti_forensics_a0_a1_a2_demo(authentic_wav_path, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load and preprocess authentic trace.
    fprintf('\nINFO[AF]: Loading and resampling authentic trace.\n');
    do_unit_peak_normalize = true;
    [x_authentic, fs_authentic] = af_load_mono_trace_resampled(authentic_wav_path, cfg.ANALYSIS_FS_HZ, do_unit_peak_normalize);

    if isempty(x_authentic)
        error('ERROR: Failed to load authentic trace for anti-forensics demo.');
    end

    fs = fs_authentic;

    fprintf('INFO[AF]: Loaded and resampled authentic trace to %d Hz. Duration = %.2f s\n', ...
            fs, numel(x_authentic)/fs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    runA0_2 = af_cfg_has_attack_technique(cfg, 'A0_2 ENF removal narrow cascade');
    runA1_2 = af_cfg_has_attack_technique(cfg, 'A1_2 ENF noise fill narrow cascade');
    runA2_2 = af_cfg_has_attack_technique(cfg, 'A2_2 ENF embedding narrow cascade');

    if ~(runA0_2 || runA1_2 || runA2_2)
        warning('AF:noAttacksEnabled', ['AF_CFG.ATTACK_TECHNIQUES contains no enabled attacks. ' ...
                                        'Nothing to do.']);
        return;
    end

    fprintf(['\nINFO[AF]: Enabled Attack techniques:' ...
             '\nA0_2 ENF removal narrow cascade:%d' ...
             '\nA1_2 ENF noise fill narrow cascade:%d' ...
             '\nA2_2 ENF embedding narrow cascade:%d\n'], ...
            runA0_2, runA1_2, runA2_2);

    runD0 = af_cfg_has_detection_technique(cfg, 'D0 enf correlation check');
    runD1 = af_cfg_has_detection_technique(cfg, 'D1 inter harmonic consistency check');
    fprintf(['\nINFO[AF]: Enabled Detection techniques:  \nD0 enf correlation check:%d ' ...
             '\nD1 inter harmonic consistency check:%d\n'], ...
            runD0, runD1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize variables for attack results (forged traces) and intermediate data.
    x_case_A0_2 = [];
    x_case_A1_2 = [];
    x_case_A2_2 = [];
    n_cascade_passes = max(1, round(af_get_cfg_numeric(cfg, 'A0_CASCADE_PASSES', 2)));

    fprintf('\nINFO[AF]: Generating forged traces for enabled attack techniques.\n');

    % ---------------------------------------------------------------
    % A0_2 - Narrow band, cascade (deeper attenuation, same spectral footprint)
    % ---------------------------------------------------------------
    if runA0_2
        fprintf('\nINFO[AF][A0_2 ENF Removal narrow cascade (%d passes)]: START\n', n_cascade_passes);
        af_print_a0_bandstop_filter_design_spec(fs, cfg);
        x_tmp = x_authentic;
        for iPass = 1:n_cascade_passes
            x_tmp = af_apply_enf_bandstop_removal(x_tmp, fs, cfg);
        end
        x_case_A0_2 = x_tmp;
        fprintf('INFO[AF][A0_2 ENF Removal narrow cascade]: END\n');
    end

    % ---------------------------------------------------------------
    % Ensure exact A0_2 base exists for A1_2/A2_2 variants.
    % ---------------------------------------------------------------
    use_exact_a0_variant_base = af_get_cfg_logical(cfg, 'A1_USE_EXACT_A0_VARIANT_BASE', true);
    if use_exact_a0_variant_base
        if (runA1_2 || runA2_2) && isempty(x_case_A0_2)
            x_tmp = x_authentic;
            for iPass = 1:n_cascade_passes
                x_tmp = af_apply_enf_bandstop_removal(x_tmp, fs, cfg);
            end
            x_case_A0_2 = x_tmp;
        end
    end

    % ---------------------------------------------------------------
    % A1 family - noise fill matched to A0_2
    %   A1_2 uses the A0_2 base plus narrow noise fill.
    % ---------------------------------------------------------------
    x_noise_fill_narrow = [];

    if runA1_2 || runA2_2
        fprintf('INFO[AF]: Computing narrow-band noise fill for A1_2/A2_2.\n');
        [~, x_noise_fill_narrow] = af_build_a1_from_a0_variant(x_authentic, x_case_A0_2, fs, cfg);
    end

    if runA1_2
        fprintf('\nINFO[AF][A1_2 ENF Noise Fill narrow cascade]: START\n');
        x_case_A1_2 = x_case_A0_2 + x_noise_fill_narrow;
        fprintf('INFO[AF][A1_2 ENF Noise Fill narrow cascade]: noise fill RMS = %.6g\n', ...
            sqrt(mean(x_noise_fill_narrow.^2, 'omitnan')));
        fprintf('INFO[AF][A1_2 ENF Noise Fill narrow cascade]: END\n');
    end

    % ---------------------------------------------------------------
    % A2_2 - A1_2 base + peak-matched FM syn ENF (narrow cascade)
    %
    % Magnitude M chosen by binary search: (A1_2 + M*syn) has same FFT peak
    % at nominal ENF freq as the authentic trace.
    % ---------------------------------------------------------------

    % Generate narrow syn ENF unit signal for A2_2
    x_syn_unit_narrow = [];
    if runA2_2
        fprintf('INFO[AF]: Generating narrow-band FM syn ENF unit signal for A2_2.\n');
        x_syn_unit_narrow = af_generate_syn_enf_component(numel(x_authentic), fs, cfg);
    end

    % Build A1_2 base needed for A2_2 (if not already computed by A1 family above)
    if runA2_2 && isempty(x_case_A1_2)
        x_case_A1_2 = x_case_A0_2 + x_noise_fill_narrow;
    end

    if runA2_2
        fprintf('\nINFO[AF][A2_2 ENF Embedding narrow cascade]: START\n');
        [M_a2_2, mi_a2_2] = af_a1_match_syn_peak_magnitude(x_authentic, x_case_A1_2, x_syn_unit_narrow, fs, cfg);
        x_case_A2_2 = x_case_A1_2 + M_a2_2 * x_syn_unit_narrow;
        fprintf(['INFO[AF][A2_2 ENF Embedding narrow cascade]: M=%.6g, target peak=%.4g, ' ...
                 'base peak=%.4g, final peak=%.4g (at %.2f Hz +/- %.2f Hz)\n'], ...
                M_a2_2, mi_a2_2.target_peak_abs, mi_a2_2.base_peak_abs, mi_a2_2.final_peak_abs, ...
                mi_a2_2.peak_match_freq_hz, mi_a2_2.peak_match_bw_hz);
        fprintf('INFO[AF][A2_2 ENF Embedding narrow cascade]: END\n');
    end

    fprintf('\nINFO[AF]: Generating consolidated FFT magnitude/attenuation comparison report.\n');
    af_print_fft_harmonic_comparison_report_combined(x_authentic, ...
        x_case_A0_2, x_case_A1_2, x_case_A2_2, ...
        fs, cfg, runA0_2, runA1_2, runA2_2);

    fprintf('\nINFO[AF]: Generating attack plots.\n');
    af_plot_attack_results(x_authentic, x_case_A0_2, x_case_A1_2, x_case_A2_2, ...
                           fs, cfg, runA0_2, runA1_2, runA2_2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nINFO[AF]: Running enabled detection techniques on the forged traces.\n');

    % D0: ENF correlation check (temporally aligned ENF extraction + overlay plotting).
    if runD0
        af_log_detector_start('Running D0 ENF correlation check on the forged traces.', 'D0 ENF correlation check', 'START');
        af_enf_results = struct();
        af_enf_results.analysis_fs_hz = fs;
        af_enf_results.enf_freq_hz = af_get_cfg_numeric(cfg, 'D0_ENF_FREQ_HZ', cfg.NOMINAL_FREQ_HZ);
        af_enf_results.frame_size_samples = af_get_cfg_numeric(cfg, 'D0_FRAME_SIZE_SAMPLES', 8000);
        af_enf_results.overlap_size_samples = af_get_cfg_numeric(cfg, 'D0_OVERLAP_SAMPLES', 0);

        fprintf('INFO[AF][D0 ENF correlation check]: Extracting ENF signature and generating comparison plot for each enabled attack technique.\n');
        d0_specs = af_get_d0_pair_specs(runA0_2, runA1_2, runA2_2);
        for iPair = 1:numel(d0_specs)
            spec = d0_specs(iPair);
            x_ref = x_authentic;
            x_forged_pair = af_get_forged_signal_by_attack_id(spec.attack_id, x_case_A0_2, x_case_A1_2, x_case_A2_2);
            pair_out = af_extract_method1_enf_pair( ...
                spec.case_label, x_ref, x_forged_pair, fs, cfg, ...
                'D0', spec.attack_id, spec.comparison_id);
            af_enf_results = af_set_struct_field_by_path(af_enf_results, spec.result_key, pair_out);
        end

        af_log_detector_end('D0 enf correlation check', 'END');
    end

    % Inter-frequency consistency check: compare ENF extracted at multiple
    % harmonics from the same signal (e.g., 60 Hz vs 120 Hz) after normalization.
    if runD1
        af_log_detector_start('Running D1 inter harmonic consistency check on the forged traces.', 'D1 inter harmonic consistency check', 'START.');
        fprintf('INFO[AF][D1 inter harmonic consistency check]: Comparing ENF phase/frequency extracted at multiple harmonics for each enabled attack technique.\n');
        af_run_harmonic_consistency_detection( ...
            x_authentic, x_case_A0_2, x_case_A1_2, x_case_A2_2, ...
            fs, cfg, runA0_2, runA1_2, runA2_2);
        af_log_detector_end('D1 inter harmonic consistency check', 'END');
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [log_file_path, plot_export_dir, prev_fig_visibility] = af_prepare_run_logging_and_plot_export(cfg, script_name, script_dir)
    if nargin < 2 || strlength(string(script_name)) == 0
        script_name = "enf_analysis_top_pair_wise_af";
    end
    if nargin < 3 || isempty(script_dir)
        script_dir = pwd;
    end
    script_name = string(script_name);
    script_dir = char(script_dir);
    log_dir_raw = char(af_cfg_get_string_or_default(cfg, 'LOG_DIR', fullfile('exp_logs')));
    plot_export_dir_raw = char(af_cfg_get_string_or_default(cfg, 'PLOT_EXPORT_DIR', fullfile('exp_results', 'af_analysis_pairwise')));
    log_dir = af_resolve_dir_path(log_dir_raw, script_dir);
    plot_export_dir = af_resolve_dir_path(plot_export_dir_raw, script_dir);

    if ~exist(log_dir, 'dir')
        [ok_log, msg_log] = mkdir(log_dir);
        if ~ok_log && ~exist(log_dir, 'dir')
            error('ERROR: Failed to create log directory: %s (%s)', log_dir, msg_log);
        end
    end
    if ~exist(plot_export_dir, 'dir')
        [ok_plot, msg_plot] = mkdir(plot_export_dir);
        if ~ok_plot && ~exist(plot_export_dir, 'dir')
            error('ERROR: Failed to create plot export directory: %s (%s)', plot_export_dir, msg_plot);
        end
    end

    log_file_path = fullfile(log_dir, sprintf('%s_log', char(script_name)));
    try
        diary off;
    catch
        % no-op
    end
    if exist(log_file_path, 'file')
        try
            delete(log_file_path);
        catch ME
            warning('AF:LogFileDeleteFailed', ...
                    'Could not delete old log file "%s": %s', log_file_path, ME.message);
        end
    end
    diary(log_file_path);

    prev_fig_visibility = get(groot, 'DefaultFigureVisible');
    if ~af_get_cfg_logical(cfg, 'PLOT_DISPLAY_ENABLE', false)
        set(groot, 'DefaultFigureVisible', 'off');
    end
end


function af_finalize_run_logging_and_plot_export(prev_fig_visibility)
    if nargin >= 1 && ~isempty(prev_fig_visibility)
        try
            set(groot, 'DefaultFigureVisible', prev_fig_visibility);
        catch
            % no-op
        end
    end
    try
        diary off;
    catch
        % no-op
    end
end


function af_log_detector_start(info_line, tag_line, start_token)
    if nargin < 3 || isempty(start_token)
        start_token = 'START';
    end
    fprintf('\nINFO[AF]: %s\n', char(string(info_line)));
    fprintf('INFO[AF][%s]: %s\n', char(string(tag_line)), char(string(start_token)));
end

function af_log_detector_end(tag_line, end_token)
    if nargin < 2 || isempty(end_token)
        end_token = 'END';
    end
    fprintf('INFO[AF][%s]: %s\n', char(string(tag_line)), char(string(end_token)));
end


function af_save_figure_svg_and_close(fh, cfg, out_name_override)
    if isempty(fh) || ~ishandle(fh)
        return;
    end
    if nargin < 3
        out_name_override = '';
    end

    plot_export_dir = af_cfg_get_string_or_default(cfg, 'RUNTIME_PLOT_EXPORT_DIR', ...
        af_cfg_get_string_or_default(cfg, 'PLOT_EXPORT_DIR', fullfile('exp_results', 'af_analysis_pairwise')));
    plot_export_dir = char(plot_export_dir);
    if ~exist(plot_export_dir, 'dir')
        [ok_plot, msg_plot] = mkdir(plot_export_dir);
        if ~ok_plot && ~exist(plot_export_dir, 'dir')
            warning('AF:PlotExportDirCreateFailed', ...
                    'Failed to create plot export directory: %s (%s).', plot_export_dir, msg_plot);
            try
                close(fh);
            catch
                % no-op
            end
            return;
        end
    end

    if ~isempty(out_name_override)
        [~, base_name, ext_name] = fileparts(char(string(out_name_override)));
        if isempty(base_name)
            base_name = 'plot';
        end
        base_name = af_sanitize_filename(base_name);
        if isempty(ext_name) || ~strcmpi(ext_name, '.png')
            ext_name = '.png';
        end
        out_name = [base_name, ext_name];
    else
        [prefix_code, short_slug] = af_build_plot_export_name_parts(fh);
        prefix_code = af_sanitize_filename(prefix_code);
        if isempty(prefix_code)
            prefix_code = 'PLOT';
        end

        short_slug = af_sanitize_filename(short_slug);
        if isempty(short_slug)
            short_slug = 'plot';
        end

        max_full_path_len = 220; % conservative for Windows MAX_PATH-sensitive environments
        max_slug_len = max(12, max_full_path_len - numel(plot_export_dir) - 1 - numel(prefix_code) - 1 - numel('.png'));
        if numel(short_slug) > max_slug_len
            short_slug = short_slug(1:max_slug_len);
            short_slug = af_sanitize_filename(short_slug);
            if isempty(short_slug)
                short_slug = 'plot';
            end
        end
        out_name = sprintf('%s_%s.png', prefix_code, short_slug);
    end
    out_path = fullfile(plot_export_dir, out_name);

    try
        if exist(out_path, 'file')
            delete(out_path);
        end
        saveas(fh, out_path);
        fprintf('INFO[AF][PlotExport]: Saved %s\n', out_path);
    catch ME
        warning('AF:PlotSaveFailed', 'Failed to save plot "%s": %s', out_path, ME.message);
    end

    svg_path = regexprep(out_path, '\.png$', '.svg', 'ignorecase');
    try
        if exist(svg_path, 'file')
            delete(svg_path);
        end
        saveas(fh, svg_path, 'svg');
        fprintf('INFO[AF][PlotExport]: Saved %s\n', svg_path);
    catch ME
        warning('AF:PlotSaveFailed', 'Failed to save plot "%s": %s', svg_path, ME.message);
    end

    try
        close(fh);
    catch
        % no-op
    end
end


function [prefix_code, short_slug] = af_build_plot_export_name_parts(fh)
    fig_name = string(get(fh, 'Name'));
    fig_title = af_get_figure_primary_title(fh);
    all_txt = lower(char(strtrim(strjoin([fig_name, fig_title], " "))));
    prefix_code = af_detect_plot_prefix_code(all_txt);

    src = strtrim(fig_name);
    if strlength(src) == 0
        src = strtrim(fig_title);
    end
    if strlength(src) == 0
        src = "plot";
    end

    parts = split(src, '|');
    if numel(parts) >= 2
        rhs = strtrim(parts(end));
        if strlength(rhs) >= 4
            src = rhs;
        else
            src = strtrim(parts(end-1) + "_" + parts(end));
        end
    end

    src = regexprep(char(src), '(?i)^\s*(A[0-2]|D[0-4])[_\-\s:]+', '');
    src = regexprep(src, '(?i)\benf_analysis_top_pair_wise_af\b', '');
    src = strtrim(string(src));
    if strlength(src) == 0
        src = "plot";
    end
    short_slug = lower(char(src));
end


function prefix_code = af_detect_plot_prefix_code(all_txt_lower)
    t = lower(char(all_txt_lower));
    prefix_code = 'PLOT';

    % Prefer detector IDs when both detector and attack hints are present.
    if contains(t, 'd1') || contains(t, 'inter harmonic') || contains(t, 'harmonic consistency') || ...
       contains(t, 'phase consistency') || contains(t, 'phase compare') || contains(t, 'phasecompare')
        prefix_code = 'D1'; return;
    end
    if contains(t, 'd0') || contains(t, 'enf correlation') || contains(t, 'temporally aligned enf')
        prefix_code = 'D0'; return;
    end

    if contains(t, 'a0') || contains(t, 'enf removal') || contains(t, 'bandstop')
        prefix_code = 'A0'; return;
    end
    if contains(t, 'a1') || contains(t, 'noise fill')
        prefix_code = 'A1'; return;
    end
    if contains(t, 'a2') || contains(t, 'embedding synthetic') || contains(t, 'syn enf')
        prefix_code = 'A2'; return;
    end
end


function title_txt = af_get_figure_primary_title(fh)
    title_txt = "";
    axes_list = findall(fh, 'Type', 'axes');
    if isempty(axes_list)
        return;
    end
    for iA = 1:numel(axes_list)
        try
            t = get(get(axes_list(iA), 'Title'), 'String');
        catch
            continue;
        end
        s = strtrim(af_any_text_to_string(t));
        if strlength(s) > 0
            title_txt = s;
            return;
        end
    end
end


function s = af_any_text_to_string(v)
    if isstring(v)
        s = strjoin(v(:).', " ");
        return;
    end
    if ischar(v)
        s = string(v);
        return;
    end
    if iscell(v)
        ss = strings(1, numel(v));
        for iV = 1:numel(v)
            ss(iV) = af_any_text_to_string(v{iV});
        end
        s = strjoin(ss, " ");
        return;
    end
    s = "";
end


function af_print_table_records_as_columns(T, record_label, col_headers, metric_var_names)
    if nargin < 2 || isempty(record_label)
        record_label = 'Record';
    end
    if nargin < 3
        col_headers = {};
    end
    if nargin < 4
        metric_var_names = {};
    end
    if isempty(T) || height(T) == 0
        fprintf('  No records.\n');
        return;
    end

    var_names = T.Properties.VariableNames;
    if isempty(metric_var_names)
        metric_var_names = var_names;
    else
        metric_var_names = cellstr(string(metric_var_names));
        metric_var_names = metric_var_names(ismember(metric_var_names, var_names));
    end
    if isempty(metric_var_names)
        fprintf('  No metrics to print.\n');
        return;
    end

    n_records = height(T);
    n_metrics = numel(metric_var_names);
    rec_headers = cell(1, n_records);
    if ~isempty(col_headers)
        tmp = cellstr(string(col_headers(:).'));
        if numel(tmp) == n_records
            rec_headers = tmp;
        end
    end
    if any(cellfun(@isempty, rec_headers))
        for j = 1:n_records
            rec_headers{j} = sprintf('%s%d', record_label, j);
        end
    end

    metric_col_w = max(numel('Metric'), max(cellfun(@numel, metric_var_names)));
    rec_col_w = zeros(1, n_records);
    for j = 1:n_records
        rec_col_w(j) = numel(rec_headers{j});
    end

    values = cell(n_metrics, n_records);
    for i = 1:n_metrics
        vn = metric_var_names{i};
        col = T.(vn);
        for j = 1:n_records
            s = af_table_col_value_to_string(col, j);
            values{i, j} = s;
            rec_col_w(j) = max(rec_col_w(j), numel(s));
        end
    end

    fprintf('  Total records: %d\n', n_records);
    header = sprintf('%-*s', metric_col_w, 'Metric');
    for j = 1:n_records
        header = [header, ' | ', sprintf('%-*s', rec_col_w(j), rec_headers{j})]; %#ok<AGROW>
    end
    fprintf('  %s\n', header);

    sep_len = metric_col_w + sum(rec_col_w) + (3 * n_records);
    fprintf('  %s\n', repmat('-', 1, sep_len));

    for i = 1:n_metrics
        line = sprintf('%-*s', metric_col_w, metric_var_names{i});
        for j = 1:n_records
            line = [line, ' | ', sprintf('%-*s', rec_col_w(j), values{i, j})]; %#ok<AGROW>
        end
        fprintf('  %s\n', line);
    end
end


function s = af_table_col_value_to_string(col, iRow)
    if iscell(col)
        try
            v = col{iRow};
        catch
            v = col(iRow);
        end
    else
        try
            v = col(iRow, :);
        catch
            try
                v = col(iRow);
            catch
                v = col;
            end
        end
    end

    if isstring(v)
        if isscalar(v)
            s = char(v);
        else
            s = strjoin(cellstr(v(:).'), ' | ');
        end
        s = strtrim(s);
        return;
    end
    if ischar(v)
        s = strtrim(v);
        return;
    end
    if iscategorical(v)
        s = char(string(v));
        s = strtrim(s);
        return;
    end
    if islogical(v)
        if isscalar(v)
            s = char(string(v));
        else
            parts = arrayfun(@(x) char(string(x)), v(:).', 'UniformOutput', false);
            s = ['[', strjoin(parts, ' '), ']'];
        end
        return;
    end
    if isnumeric(v)
        vv = double(v);
        if isempty(vv)
            s = '[]';
            return;
        end
        if isscalar(vv)
            s = sprintf('%.6g', vv);
            return;
        end
        parts = arrayfun(@(x) sprintf('%.6g', x), vv(:).', 'UniformOutput', false);
        s = ['[', strjoin(parts, ' '), ']'];
        return;
    end
    try
        s = char(string(v));
        s = strtrim(regexprep(s, '\s+', ' '));
    catch
        s = '<unprintable>';
    end
end


function out_path = af_resolve_dir_path(path_in, base_dir)
    path_in = char(string(path_in));
    base_dir = char(string(base_dir));
    if af_is_absolute_path(path_in)
        out_path = path_in;
    else
        out_path = fullfile(base_dir, path_in);
    end
end


function tf = af_is_absolute_path(path_in)
    path_in = char(string(path_in));
    if isempty(path_in)
        tf = false;
        return;
    end
    if ispc
        tf = ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]|^\\\\', 'once'));
    else
        tf = startsWith(path_in, '/');
    end
end


function out_name = af_sanitize_filename(in_name)
    out_name = char(string(in_name));
    out_name = regexprep(out_name, '[^A-Za-z0-9_-]+', '_');
    out_name = regexprep(out_name, '_+', '_');
    out_name = regexprep(out_name, '^_+', '');
    out_name = regexprep(out_name, '_+$', '');
    if isempty(out_name)
        out_name = 'figure';
    end
end



function tf = af_get_cfg_logical(cfg, field_name, default_val)
    tf = default_val;
    if isfield(cfg, field_name)
        try
            tf = logical(cfg.(field_name));
        catch
            tf = default_val;
        end
    end
    tf = logical(tf);
end


function val = af_get_cfg_numeric(cfg, field_name, default_val)
    val = default_val;
    if isfield(cfg, field_name)
        try
            v = double(cfg.(field_name));
            if isscalar(v) && isfinite(v)
                val = v;
            end
        catch
            val = default_val;
        end
    end
end


function tf = af_cfg_has_attack_technique(cfg, attack_name)
    tf = false;
    if nargin < 2 || isempty(attack_name)
        return;
    end
    want = af_cfg_canonical_attack_name(attack_name);

    items = {};
    if isfield(cfg, 'ATTACK_TECHNIQUES') && ~isempty(cfg.ATTACK_TECHNIQUES)
        raw = cfg.ATTACK_TECHNIQUES;
        try
            if ischar(raw)
                items = {raw};
            elseif isstring(raw)
                items = cellstr(raw(:));
            elseif iscell(raw)
                items = raw(:).';
            end
        catch
            items = {};
        end
    else
        % Default: enabled attack techniques (_2 variants only).
        items = {'A0_2 ENF removal narrow cascade', ...
                 'A1_2 ENF noise fill narrow cascade', ...
                 'A2_2 ENF embedding narrow cascade'};
    end

    for i = 1:numel(items)
        try
            s = af_cfg_canonical_attack_name(items{i});
            if strcmp(s, want)
                tf = true;
                return;
            end
        catch
            % ignore malformed entry
        end
    end
end


function s = af_cfg_canonical_attack_name(name_in)
    s = lower(strtrim(char(name_in)));
    s = strrep(s, '-', '_');
    s = strrep(s, ' ', '_');
    % A0_2 canonical name
    if strcmp(s, 'a0_2') || strcmp(s, 'a0_2_enf_removal_narrow_cascade')
        s = 'a0_2_enf_removal_narrow_cascade';
    % A1_2 canonical name
    elseif strcmp(s, 'a1_2') || strcmp(s, 'a1_2_enf_noise_fill_narrow_cascade')
        s = 'a1_2_enf_noise_fill_narrow_cascade';
    % A2_2 canonical name
    elseif strcmp(s, 'a2_2') || strcmp(s, 'a2_2_enf_embedding_narrow_cascade')
        s = 'a2_2_enf_embedding_narrow_cascade';
    end
end


function tf = af_cfg_has_detection_technique(cfg, technique_name)
    % Detection-technique selector driven by AF_CFG.DETECTION_TECHNIQUES.
    % Supports current D0..D1 names and legacy aliases.
    tf = false;
    if nargin < 2 || isempty(technique_name)
        return;
    end
    want = af_cfg_canonical_detection_name(technique_name);
    items = {};
    if isfield(cfg, 'DETECTION_TECHNIQUES') && ~isempty(cfg.DETECTION_TECHNIQUES)
        raw = cfg.DETECTION_TECHNIQUES;
        try
            if ischar(raw)
                items = {raw};
            elseif isstring(raw)
                items = cellstr(raw(:));
            elseif iscell(raw)
                items = raw(:).';
            end
        catch
            items = {};
        end
    else
        % Default: enabled detectors.
        items = { ...
            'D0 enf correlation check', ...
            'D1 inter harmonic consistency check' ...
        };
    end

    for i = 1:numel(items)
        try
            s = af_cfg_canonical_detection_name(items{i});
            if strcmp(s, want)
                tf = true;
                return;
            end
        catch
            % ignore malformed entry
        end
    end
end


function s = af_cfg_canonical_detection_name(name_in)
    s = lower(strtrim(char(name_in)));
    s = strrep(s, '-', '_');
    s = strrep(s, ' ', '_');
    if strcmp(s, 'harmonic_consistency')
        s = 'd1_inter_harmonic_consistency_check';
    elseif strcmp(s, 'd0')
        s = 'd0_enf_correlation_check';
    elseif strcmp(s, 'd1')
        s = 'd1_inter_harmonic_consistency_check';
    end
end


function attack_ids = af_get_enabled_attack_ids(runA0_2, runA1_2, runA2_2)
    if nargin < 1, runA0_2 = false; end
    if nargin < 2, runA1_2 = false; end
    if nargin < 3, runA2_2 = false; end
    attack_ids = {};
    if runA0_2, attack_ids{end+1} = 'A0_2'; end
    if runA1_2, attack_ids{end+1} = 'A1_2'; end
    if runA2_2, attack_ids{end+1} = 'A2_2'; end
end


function meta = af_get_attack_meta(attack_id)
    aid = upper(strtrim(char(string(attack_id))));
    switch aid
        case 'A0_2'
            meta = struct( ...
                'id', 'A0_2', ...
                'plot_case_label', 'A0_2 ENF Removal (Narrow Bandstop, Cascade)', ...
                'd0_case_label_authentic', 'A0_2 ENF Removal Narrow Cascade - Authentic vs Forged', ...
                'd1_cfg_key', 'D1_PLOT_A0_2', 'd1_cfg_default', false, ...
                'd1_key', 'A0_2_enf_removal_narrow_cascade_forged_trace', ...
                'd1_label', 'A0_2 - ENF Removal Narrow Cascade (Forged trace)', ...
                'd1_trace_label', 'A0_2 Forged trace');
        case 'A1_2'
            meta = struct( ...
                'id', 'A1_2', ...
                'plot_case_label', 'A1_2 ENF Noise Fill (Narrow Bandstop, Cascade)', ...
                'd0_case_label_authentic', 'A1_2 ENF Noise Fill Narrow Cascade - Authentic vs Forged', ...
                'd1_cfg_key', 'D1_PLOT_A1', 'd1_cfg_default', false, ...
                'd1_key', 'A1_2_enf_noise_fill_narrow_cascade_forged_trace', ...
                'd1_label', 'A1_2 - ENF Noise Fill Narrow Cascade (Forged trace)', ...
                'd1_trace_label', 'A1_2 Forged trace');
        case 'A2_2'
            meta = struct( ...
                'id', 'A2_2', ...
                'plot_case_label', 'A2_2 ENF Synthetic Embedding (Narrow, Cascade)', ...
                'd0_case_label_authentic', 'A2_2 ENF Synthetic Embedding Narrow Cascade - Authentic vs Forged', ...
                'd1_cfg_key', 'D1_PLOT_A2', 'd1_cfg_default', true, ...
                'd1_key', 'A2_2_enf_embedding_narrow_cascade_forged_trace', ...
                'd1_label', 'A2_2 - Synthetic ENF Embedding Narrow Cascade (Forged trace)', ...
                'd1_trace_label', 'A2_2 Forged trace');
        otherwise
            error('AF:UnknownAttackID', 'Unknown attack id: %s', aid);
    end
end


% af_get_d0_pair_specs - Return D0 authentic-vs-forged pair specifications for enabled attacks.
%
%   Convenience wrapper around af_get_detector_pair_specs for detector 'D0'.
%
%   Inputs:
%     runA0_2 - Logical; include A0_2 pair if true.
%     runA1_2 - Logical; include A1_2 pair if true.
%     runA2_2 - Logical; include A2_2 pair if true.
%
%   Outputs:
%     specs - Struct array with fields: result_key, case_label, attack_id,
%             comparison_id, label_unforged, label_forged.
%
function specs = af_get_d0_pair_specs(runA0_2, runA1_2, runA2_2)
    if nargin < 1, runA0_2 = false; end
    if nargin < 2, runA1_2 = false; end
    if nargin < 3, runA2_2 = false; end
    specs = af_get_detector_pair_specs('D0', runA0_2, runA1_2, runA2_2);
end



% af_get_detector_pair_specs - Build pair specification structs for a given detector.
%
%   Iterates over enabled attack IDs and constructs one spec struct per pair,
%   setting the result storage key, case label, attack ID, and comparison ID
%   appropriate for the specified detector (currently supports 'D0').
%
%   Inputs:
%     detector_id - Detector identifier string (e.g., 'D0').
%     runA0_2     - Logical; include A0_2 if true.
%     runA1_2     - Logical; include A1_2 if true.
%     runA2_2     - Logical; include A2_2 if true.
%
%   Outputs:
%     specs - Struct array with one entry per enabled attack x detector pair.
%
function specs = af_get_detector_pair_specs(detector_id, runA0_2, runA1_2, runA2_2)
    if nargin < 2, runA0_2 = false; end
    if nargin < 3, runA1_2 = false; end
    if nargin < 4, runA2_2 = false; end
    did = upper(strtrim(char(string(detector_id))));
    attack_ids = af_get_enabled_attack_ids(runA0_2, runA1_2, runA2_2);
    empty_spec = struct('result_key', '', 'case_label', '', 'attack_id', '', ...
                        'comparison_id', '', 'label_unforged', '', 'label_forged', '');
    specs = repmat(empty_spec, 1, numel(attack_ids));
    spec_count = 0;
    for i = 1:numel(attack_ids)
        meta = af_get_attack_meta(attack_ids{i});
        switch did
            case 'D0'
                spec_count = spec_count + 1;
                specs(spec_count) = struct( ...
                    'result_key', sprintf('attack%s', meta.id), ...
                    'case_label', meta.d0_case_label_authentic, ...
                    'attack_id', meta.id, ...
                    'comparison_id', 'authentic_vs_forged', ...
                    'label_unforged', '', ...
                    'label_forged', '');
            otherwise
                error('AF:UnknownDetectorID', 'Unknown detector id for pair specs: %s', did);
        end
    end
    specs = specs(1:spec_count);
end


% af_get_forged_signal_by_attack_id - Return the forged signal for a given attack ID.
%
%   Inputs:
%     attack_id   - Attack ID string: 'A0_2', 'A1_2', or 'A2_2'.
%     x_case_A0_2 - Forged signal from A0_2 (or [] if not computed).
%     x_case_A1_2 - Forged signal from A1_2 (or [] if not computed).
%     x_case_A2_2 - Forged signal from A2_2 (or [] if not computed).
%
%   Outputs:
%     x - Selected forged signal vector, or [] for an unknown attack ID.
%
function x = af_get_forged_signal_by_attack_id(attack_id, x_case_A0_2, x_case_A1_2, x_case_A2_2)
    if nargin < 2, x_case_A0_2 = []; end
    if nargin < 3, x_case_A1_2 = []; end
    if nargin < 4, x_case_A2_2 = []; end
    aid = upper(strtrim(char(string(attack_id))));
    switch aid
        case 'A0_2', x = x_case_A0_2;
        case 'A1_2', x = x_case_A1_2;
        case 'A2_2', x = x_case_A2_2;
        otherwise,   x = [];
    end
end


% af_set_struct_field_by_path - Set a nested struct field using a dot-separated path.
%
%   Supports one or two levels of nesting (e.g., 'top' or 'top.leaf').
%
%   Inputs:
%     s         - Input struct to modify.
%     path_expr - Dot-separated field path string (e.g., 'attackA0_2' or 'a.b').
%     value     - Value to assign to the specified field.
%
%   Outputs:
%     s - Modified struct with the specified field set to value.
%
function s = af_set_struct_field_by_path(s, path_expr, value)
    parts = regexp(char(path_expr), '\.', 'split');
    if isempty(parts)
        return;
    end
    if numel(parts) == 1
        s.(parts{1}) = value;
        return;
    end

    root = parts{1};
    leaf = parts{2};
    if ~isfield(s, root) || ~isstruct(s.(root)) || isempty(s.(root))
        s.(root) = struct();
    end
    tmp = s.(root);
    tmp.(leaf) = value;
    s.(root) = tmp;
end


% af_extract_method1_enf_pair - Extract ENF from an authentic/forged pair and report D0 results.
%
%   Calls proc_enf_analysis with lag optimisation disabled to extract ENF
%   signatures from both the authentic and forged signals, computes Pearson
%   correlation, saves an overlay plot, and returns a result struct.
%
%   Inputs:
%     case_label    - Descriptive label for this pair (used in plot title and log).
%     x_authentic   - Authentic signal vector (samples at fs Hz).
%     x_forged      - Forged signal vector (samples at fs Hz).
%     fs            - Sampling rate in Hz.
%     cfg           - AF_CFG struct with D0_FRAME_SIZE_SAMPLES, D0_OVERLAP_SAMPLES,
%                     D0_ENF_FREQ_HZ, RUNTIME_PLOT_EXPORT_DIR, NOMINAL_FREQ_HZ.
%     detector_id   - (optional) Detector tag string for file naming (default 'D0').
%     attack_id     - (optional) Attack tag string for file naming (default 'A0').
%     comparison_id - (optional) Comparison tag for file naming (default 'authentic_vs_forged').
%
%   Outputs:
%     out - Struct with fields: case_label, method, enf_freq_hz, frame_size_samples,
%           overlap_size_samples, fs_hz, hop_sec, corr0, matched_corr, bestLagSec,
%           bestLagFrames, time_authentic_trace_sec, time_forged_trace_sec,
%           enf_authentic_trace_hz, enf_forged_trace_hz.
%
function out = af_extract_method1_enf_pair(case_label, x_authentic, x_forged, fs, cfg, detector_id, attack_id, comparison_id)
    if nargin < 6 || isempty(detector_id)
        detector_id = 'D0';
    end
    if nargin < 7 || isempty(attack_id)
        attack_id = 'A0';
    end
    if nargin < 8 || isempty(comparison_id)
        comparison_id = 'authentic_vs_forged';
    end
    frame_size = round(af_get_cfg_numeric(cfg, 'D0_FRAME_SIZE_SAMPLES', 8000));
    overlap_size = round(af_get_cfg_numeric(cfg, 'D0_OVERLAP_SAMPLES', 0));
    enf_freq_hz = af_get_cfg_numeric(cfg, 'D0_ENF_FREQ_HZ', cfg.NOMINAL_FREQ_HZ);
    enable_plot = true;
    plot_export_dir = af_cfg_get_string_or_default(cfg, 'RUNTIME_PLOT_EXPORT_DIR', fullfile('exp_results', 'af_analysis_pairwise'));
    plot_file_stem = sprintf('%s_%s_%s_enf_overlay.png', char(detector_id), char(attack_id), char(comparison_id));

    if frame_size <= 0
        error('AF ENF Method-1 extraction requires frame_size > 0.');
    end
    if overlap_size < 0 || overlap_size >= frame_size
        error('AF ENF Method-1 extraction requires 0 <= overlap_size < frame_size.');
    end
    if numel(x_authentic) < frame_size || numel(x_forged) < frame_size
        error('AF ENF Method-1 extraction requires traces of at least frame_size samples.');
    end

    fprintf('INFO[AF][D0_enf_correlation_check]: Extracting Method-1 ENF at %.2f Hz for %s ...\n', enf_freq_hz, case_label);

    nominal_hz = af_get_cfg_numeric(cfg, 'NOMINAL_FREQ_HZ', 60);
    if isfield(cfg, 'FFT_HARMONIC_TEXT_FREQS_HZ') && ~isempty(cfg.FFT_HARMONIC_TEXT_FREQS_HZ)
        harmonics_arr = double(cfg.FFT_HARMONIC_TEXT_FREQS_HZ(:).');
    else
        harmonics_arr = nominal_hz * (1:7);
    end
    nfft_spec  = round(af_get_cfg_numeric(cfg, 'SPEC_NFFT', 32768));
    title_authentic = sprintf('%s | Authentic trace', case_label);
    title_forged = sprintf('%s | Forged trace', case_label);

    [corr0, ~, ~, ~, enf_orig, enf_mod] = proc_enf_analysis( ...
        x_authentic, x_forged, ...
        nfft_spec, frame_size, overlap_size, ...
        harmonics_arr, nominal_hz, harmonics_arr, nominal_hz, ...
        1, enf_freq_hz, [], ...
        1, enf_freq_hz, [], ...
        title_authentic, title_forged, false, ...
        'Fs1', fs, 'Fs2', fs, 'LagOptEnable', false);

    hop_sec = (frame_size - overlap_size) / fs;

    out = struct();
    out.case_label = case_label;
    out.method = 'weighted_pmf_method1';
    out.enf_freq_hz = enf_freq_hz;
    out.frame_size_samples = frame_size;
    out.overlap_size_samples = overlap_size;
    out.fs_hz = fs;
    out.hop_sec = hop_sec;
    out.corr0 = corr0;
    out.matched_corr = corr0;   % lag optimisation disabled; matched == baseline
    out.bestLagSec = 0;
    out.bestLagFrames = 0;
    out.time_authentic_trace_sec = (0:numel(enf_orig)-1).' * hop_sec;
    out.time_forged_trace_sec = (0:numel(enf_mod)-1).'  * hop_sec;
    out.enf_authentic_trace_hz = enf_orig;
    out.enf_forged_trace_hz = enf_mod;

    n_points = min(numel(enf_orig), numel(enf_mod));
    if n_points > 0
        d = enf_mod(1:n_points) - enf_orig(1:n_points);
        fprintf(['INFO[AF][D0_enf_correlation_check]: %s -> %d ENF points (hop %.3f s), ' ...
                 'Authentic mean=%.4f Hz, Forged mean=%.4f Hz, Delta RMS=%.4f Hz, Corr=%.4f\n'], ...
                case_label, n_points, hop_sec, mean(enf_orig(1:n_points)), ...
                mean(enf_mod(1:n_points)), sqrt(mean(d.^2)), corr0);
    else
        fprintf('INFO[AF][D0_enf_correlation_check]: %s -> No ENF points extracted.\n', case_label);
    end

    % Overlay plot and PNG export (plot_mode = 'overlay_only', enable_plot = true).
    if enable_plot
        f3 = figure;
        set(f3, 'Position', [890 200 540 480]);
        [t_al, s1_al, s2_al] = af_get_aligned_enf_for_plot( ...
            enf_orig, enf_mod, fs, frame_size, overlap_size, 0, 60);
        if isempty(t_al)
            warning('af_extract_method1_enf_pair:AlignedPlotSkipped', ...
                    'Aligned ENF plot skipped (insufficient overlap).');
        else
            plot(t_al, s1_al, 'b-'); grid on; hold on;
            plot(t_al, s2_al, 'r-'); hold off;
        end
        corr_str = ['(Corr = ', num2str(corr0, '%.4f'), ', Lag = 0.00 s)'];
        legend(string(title_authentic) + " ENF", string(title_forged) + " ENF", ...
               'Location', 'southoutside');
        title({'Temporally Aligned Extracted ENF Signals', corr_str});
        xlabel('Time (seconds)'); ylabel('Instantaneous ENF Value (Hz)');

        out_path = fullfile(char(plot_export_dir), char(plot_file_stem));
        out_dir  = fileparts(out_path);
        if ~isempty(out_dir) && ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        try
            saveas(f3, out_path);
            fprintf('INFO[AF][PlotExport]: Saved %s\n', out_path);
        catch ME
            warning('af_extract_method1_enf_pair:PlotSaveFailed', ...
                    'Failed to save plot "%s": %s', out_path, ME.message);
        end
        svg_path = regexprep(out_path, '\.png$', '.svg', 'ignorecase');
        try
            saveas(f3, svg_path, 'svg');
            fprintf('INFO[AF][PlotExport]: Saved %s\n', svg_path);
        catch ME
            warning('af_extract_method1_enf_pair:PlotSaveFailed', ...
                    'Failed to save plot "%s": %s', svg_path, ME.message);
        end
    end
end


% af_build_detector_scenarios - Assemble analysis scenario structs for a given detector.
%
%   For detector 'D1', builds one scenario per enabled attack plus the authentic
%   baseline, skipping forged variants whose D1 plot flag is disabled in cfg.
%
%   Inputs:
%     detector_id - Detector ID string (e.g., 'D1').
%     cfg         - AF_CFG struct with D1 plot flags and attack configuration.
%     runA0_2     - Logical; include A0_2 scenario if true.
%     runA1_2     - Logical; include A1_2 scenario if true.
%     runA2_2     - Logical; include A2_2 scenario if true.
%     x_authentic - Authentic signal vector.
%     x_case_A0_2 - A0_2 forged signal vector (or []).
%     x_case_A1_2 - A1_2 forged signal vector (or []).
%     x_case_A2_2 - A2_2 forged signal vector (or []).
%
%   Outputs:
%     scenarios - Struct array with fields: key, label, trace_label, x,
%                 x_ref_pair, ref_trace_label.
%
function scenarios = af_build_detector_scenarios(detector_id, cfg, runA0_2, runA1_2, runA2_2, x_authentic, x_case_A0_2, x_case_A1_2, x_case_A2_2)
    did = upper(strtrim(char(string(detector_id))));
    attack_ids = af_get_enabled_attack_ids(runA0_2, runA1_2, runA2_2);
    switch did
        case 'D1'
            scenario_template = struct('key', '', 'label', '', 'trace_label', '', ...
                                       'x', [], 'x_ref_pair', [], 'ref_trace_label', '');
            scenarios = repmat(scenario_template, 1, numel(attack_ids) + 1);
            scenario_count = 1;
            scenarios(scenario_count) = struct('key', 'baseline_authentic_trace', ...
                                               'label', 'Baseline (Authentic trace)', ...
                                               'trace_label', 'Authentic trace', ...
                                               'x', x_authentic, ...
                                               'x_ref_pair', [], ...
                                               'ref_trace_label', '');
            for i = 1:numel(attack_ids)
                meta = af_get_attack_meta(attack_ids{i});
                if af_get_cfg_logical(cfg, meta.d1_cfg_key, meta.d1_cfg_default)
                    scenario_count = scenario_count + 1;
                    scenarios(scenario_count) = struct( ...
                        'key', meta.d1_key, ...
                        'label', meta.d1_label, ...
                        'trace_label', meta.d1_trace_label, ...
                        'x', af_get_forged_signal_by_attack_id(meta.id, x_case_A0_2, x_case_A1_2, x_case_A2_2), ...
                        'x_ref_pair', [], ...
                        'ref_trace_label', '');
                end
            end
            scenarios = scenarios(1:scenario_count);

        otherwise
            scenarios = struct([]);
    end
end


% af_run_harmonic_consistency_detection - Run D1 inter-harmonic consistency check on all scenarios.
%
%   Extracts ENF at multiple harmonic frequencies from each scenario signal,
%   normalises each to the implied fundamental, and computes Pearson correlation
%   between the anchor harmonic and each other harmonic. Saves per-scenario plots.
%
%   Inputs:
%     x_authentic   - Authentic signal vector.
%     x_case_A0_2   - A0_2 forged signal vector (or []).
%     x_case_A1_2   - A1_2 forged signal vector (or []).
%     x_case_A2_2   - A2_2 forged signal vector (or []).
%     fs            - Sampling rate in Hz.
%     cfg           - AF_CFG struct with D1_HARMONIC_MULTS, NOMINAL_FREQ_HZ,
%                     D0_FRAME_SIZE_SAMPLES, D0_OVERLAP_SAMPLES.
%     runA0_2       - Logical; include A0_2 scenario if true.
%     runA1_2       - Logical; include A1_2 scenario if true.
%     runA2_2       - Logical; include A2_2 scenario if true.
%
%   Outputs:
%     results - Struct with fields: method, fs_hz, frame_size_samples,
%               overlap_size_samples, harmonic_mults, harmonic_pair_mults,
%               freq_pairs_hz, valid, rows, table.
%
function results = af_run_harmonic_consistency_detection(x_authentic, x_case_A0_2, x_case_A1_2, x_case_A2_2, fs, cfg, runA0_2, runA1_2, runA2_2)
    results = struct();
    results.method = 'inter_frequency_consistency_check';
    results.fs_hz = fs;
    results.frame_size_samples = round(af_get_cfg_numeric(cfg, 'D0_FRAME_SIZE_SAMPLES', 8000));
    results.overlap_size_samples = round(af_get_cfg_numeric(cfg, 'D0_OVERLAP_SAMPLES', 0));
    mults = [1 2];
    if isfield(cfg, 'D1_HARMONIC_MULTS') && ~isempty(cfg.D1_HARMONIC_MULTS)
        try
            mults = double(cfg.D1_HARMONIC_MULTS(:).');
            mults = mults(isfinite(mults) & mults > 0);
        catch
            mults = [1 2];
        end
    end
    if numel(mults) < 2
        warning('AF:HarmonicConsistencyNeedTwoHarmonics', ...
                'D1_HARMONIC_MULTS must contain at least two harmonics. Skipping.');
        results.valid = false;
        return;
    end
    % Pairing strategy: compare the first harmonic in the list (typically fundamental, e.g., 60 Hz)
    % against each additional harmonic. Example [1 2 3 4 5] -> (1,2), (1,3), (1,4), (1,5).
    % Optional sanity check: include (1,1) -> 60 vs 60.
    anchor_mult = mults(1);
    compare_mults = mults(2:end);
    compare_mults = [anchor_mult, compare_mults];
    pair_mults = [repmat(anchor_mult, numel(compare_mults), 1), compare_mults(:)];
    nominal_hz = af_get_cfg_numeric(cfg, 'NOMINAL_FREQ_HZ', 60);
    freq_pairs_hz = nominal_hz * pair_mults;
    results.harmonic_mults = mults;
    results.harmonic_pair_mults = pair_mults;
    results.freq_pairs_hz = freq_pairs_hz;
    if size(freq_pairs_hz,1) == 1
        results.freq_pair_hz = freq_pairs_hz(1,:); % backward-compatible convenience field
    else
        results.freq_pair_hz = [];
    end

    scenarios = af_build_detector_scenarios('D1', cfg, runA0_2, runA1_2, runA2_2, x_authentic, x_case_A0_2, x_case_A1_2, x_case_A2_2);

    if isempty(scenarios)
        fprintf('INFO[AF][D1_inter_harmonic_consistency_check]: No scenarios selected for harmonic consistency check.\n');
        results.valid = false;
        return;
    end

    results.valid = true;
    % Use an empty struct array without a fixed field schema because each row
    % carries additional diagnostic fields (raw ENF tracks, normalization stats, meta).
    row_results = struct([]);

    for i = 1:numel(scenarios)
        sc = scenarios(i);
        for p = 1:size(freq_pairs_hz, 1)
            fp = freq_pairs_hz(p, :);
            hp = pair_mults(p, :);
            fprintf(['INFO[AF][D1_inter_harmonic_consistency_check]: Running %g Hz vs %g Hz ENF consistency ' ...
                     '(harmonics %g vs %g) on "%s" ...\n'], ...
                    fp(1), fp(2), hp(1), hp(2), sc.label);

            rr = af_harmonic_consistency_single_signal(sc.x, fs, cfg, sc.label, sc.trace_label, fp(1), fp(2));
            rr.scenario_key = sc.key;
            rr.scenario_label = sc.label;
            rr.trace_label = sc.trace_label;
            rr.harmonic_mult_1 = hp(1);
            rr.harmonic_mult_2 = hp(2);

            if isempty(row_results)
                row_results = rr;
            else
                row_results(end+1) = rr; %#ok<AGROW>
            end

            if isfield(sc, 'x_ref_pair') && ~isempty(sc.x_ref_pair)
                rr_ref = af_harmonic_consistency_single_signal(sc.x_ref_pair, fs, cfg, sc.label, char(sc.ref_trace_label), fp(1), fp(2));
                rr_ref.scenario_key = sc.key;
                rr_ref.harmonic_mult_1 = hp(1);
                rr_ref.harmonic_mult_2 = hp(2);
                af_plot_harmonic_consistency_injection_vs_forged(rr_ref, rr, cfg);
            else
                af_plot_harmonic_consistency_overlay(rr, cfg);
            end
        end
    end

    results.rows = row_results;
    results.table = af_harmonic_consistency_make_table(row_results);

    af_harmonic_consistency_print_table(results.table);
end


% af_harmonic_consistency_single_signal - Extract ENF at two frequencies and compute consistency.
%
%   Calls proc_enf_analysis twice on the same signal (self-comparison) at
%   freq1_hz and freq2_hz, divides each series by its harmonic multiplier to
%   obtain the implied fundamental, and computes Pearson correlation between them.
%
%   Inputs:
%     x             - Signal vector (samples at fs Hz).
%     fs            - Sampling rate in Hz.
%     cfg           - AF_CFG struct with D0_FRAME_SIZE_SAMPLES, D0_OVERLAP_SAMPLES,
%                     NOMINAL_FREQ_HZ, FFT_HARMONIC_TEXT_FREQS_HZ.
%     scenario_label - Label string for this scenario (used in plot titles).
%     trace_label    - Short trace label string (used in plot legends).
%     freq1_hz       - First extraction frequency (anchor, typically fundamental).
%     freq2_hz       - Second extraction frequency (comparison harmonic).
%
%   Outputs:
%     rr - Result struct with fields: scenario_key, scenario_label, trace_label,
%          freq1_hz, freq2_hz, harmonic_mult_1/2, frame_size_samples,
%          overlap_size_samples, hop_sec, n_overlap, corr_norm, enf1_hz, enf2_hz,
%          enf1_fund_hz, enf2_fund_hz, time_sec.
%
function rr = af_harmonic_consistency_single_signal(x, fs, cfg, scenario_label, trace_label, freq1_hz, freq2_hz)
    frame_size = round(af_get_cfg_numeric(cfg, 'D0_FRAME_SIZE_SAMPLES', 8000));
    overlap_size = round(af_get_cfg_numeric(cfg, 'D0_OVERLAP_SAMPLES', 0));
    if frame_size <= 0 || overlap_size < 0 || overlap_size >= frame_size
        error('AF harmonic consistency requires valid frame/overlap settings.');
    end

    nominal_hz = af_get_cfg_numeric(cfg, 'NOMINAL_FREQ_HZ', 60);
    if isfield(cfg, 'FFT_HARMONIC_TEXT_FREQS_HZ') && ~isempty(cfg.FFT_HARMONIC_TEXT_FREQS_HZ)
        harmonics_arr = double(cfg.FFT_HARMONIC_TEXT_FREQS_HZ(:).');
    else
        harmonics_arr = nominal_hz * (1:7);
    end
    nfft_spec = round(af_get_cfg_numeric(cfg, 'SPEC_NFFT', 32768));

    [~, ~, ~, ~, enf1, enf2] = proc_enf_analysis( ...
        x, x, ...
        nfft_spec, frame_size, overlap_size, ...
        harmonics_arr, nominal_hz, harmonics_arr, nominal_hz, ...
        1, freq1_hz, [], ...
        1, freq2_hz, [], ...
        sprintf('%s | ENF @ %.1f Hz', scenario_label, freq1_hz), ...
        sprintf('%s | ENF @ %.1f Hz', scenario_label, freq2_hz), ...
        false, ...
        'Fs1', fs, 'Fs2', fs, 'LagOptEnable', false);

    hop_sec = (frame_size - overlap_size) / fs;
    n = min(numel(enf1), numel(enf2));
    enf1 = enf1(:); enf2 = enf2(:);
    t = (0:n-1).' * hop_sec;

    e1 = enf1(1:n);
    e2 = enf2(1:n);
    t = t(1:n);
    good = isfinite(e1) & isfinite(e2);
    e1 = e1(good);
    e2 = e2(good);
    t = t(good);

    % Map each harmonic measurement back to the implied fundamental frequency
    % by dividing by the harmonic number (freq_hz / nominal_hz).
    % A genuine signal satisfies: f_extracted_at_kth_harmonic / k == fundamental.
    h1 = freq1_hz / nominal_hz;
    h2 = freq2_hz / nominal_hz;
    e1_fund = e1 / h1;
    e2_fund = e2 / h2;
    corr_norm = af_calc_pearson_corr(e1_fund, e2_fund);

    rr = struct();
    rr.scenario_key = '';
    rr.scenario_label = scenario_label;
    rr.trace_label = trace_label;
    rr.freq1_hz = freq1_hz;
    rr.freq2_hz = freq2_hz;
    rr.harmonic_mult_1 = h1;
    rr.harmonic_mult_2 = h2;
    rr.frame_size_samples = frame_size;
    rr.overlap_size_samples = overlap_size;
    rr.hop_sec = (frame_size - overlap_size) / fs;
    rr.n_overlap = numel(e1_fund);
    rr.corr_norm = corr_norm;
    rr.enf1_hz = e1;
    rr.enf2_hz = e2;
    rr.enf1_fund_hz = e1_fund;
    rr.enf2_fund_hz = e2_fund;
    rr.time_sec = t;
end


% af_calc_pearson_corr - Compute Pearson correlation between two vectors.
%
%   Trims to the shorter length, removes non-finite pairs, and returns NaN
%   when fewer than two valid samples remain.
%
%   Inputs:
%     a - Numeric vector.
%     b - Numeric vector (trimmed to match length of a).
%
%   Outputs:
%     c - Pearson r scalar, or NaN if insufficient valid data.
%
function c = af_calc_pearson_corr(a, b)
    a = a(:); b = b(:);
    L = min(numel(a), numel(b));
    if L <= 1
        c = NaN;
        return;
    end
    a = a(1:L);
    b = b(1:L);
    good = isfinite(a) & isfinite(b);
    if nnz(good) <= 1
        c = NaN;
        return;
    end
    cc = corrcoef(a(good), b(good));
    if numel(cc) < 4 || ~isfinite(cc(1,2))
        c = NaN;
    else
        c = cc(1,2);
    end
end


% af_plot_harmonic_consistency_overlay - Save a two-trace implied-fundamental ENF overlay plot.
%
%   Plots enf1_fund_hz and enf2_fund_hz from a single-signal consistency result
%   struct on a shared time axis and saves the figure as PNG and SVG.
%
%   Inputs:
%     rr  - Result struct from af_harmonic_consistency_single_signal.
%     cfg - AF_CFG struct (for plot export directory resolution).
%
%   Outputs:
%     (none) - Saves D1_<scenario_slug>_<pair_slug>_enf_consistency.{png,svg}.
%
function af_plot_harmonic_consistency_overlay(rr, cfg)
    if ~isfield(rr, 'time_sec') || isempty(rr.time_sec)
        return;
    end
    c1 = [0 0.4470 0.7410];
    c2 = [0.8500 0.3250 0.0980];
    f = figure('Name', sprintf('%s - Harmonic Consistency', char(rr.scenario_label)), ...
               'Color', 'w', 'Position', [220 120 900 380]);
    ax = axes('Parent', f);
    plot(ax, rr.time_sec, rr.enf1_fund_hz, '-', 'Color', c1, 'LineWidth', 1.0); hold(ax, 'on'); grid(ax, 'on');
    plot(ax, rr.time_sec, rr.enf2_fund_hz, '--', 'Color', c2, 'LineWidth', 1.0);
    hold(ax, 'off');
    xlabel(ax, 'Time (seconds)');
    ylabel(ax, 'Implied Fundamental ENF (Hz)');
    title(ax, sprintf(['Inter-Frequency Consistency Check | %s\n' ...
                       'ENF @ %.1f Hz vs %.1f Hz (implied fundamental) | Corr = %.4f'], ...
                      char(rr.scenario_label), rr.freq1_hz, rr.freq2_hz, rr.corr_norm), ...
          'Interpreter', 'none');
    leg = legend(ax, ...
                 sprintf('Implied fundamental from %.1f Hz harmonic', rr.freq1_hz), ...
                 sprintf('Implied fundamental from %.1f Hz harmonic', rr.freq2_hz), ...
                 'Location', 'southoutside');
    try
        set(leg, 'Interpreter', 'none');
    catch
    end

    scenario_slug = 'scenario';
    if isfield(rr, 'scenario_key') && ~isempty(rr.scenario_key)
        scenario_slug = lower(af_sanitize_filename(char(rr.scenario_key)));
    end
    if isempty(scenario_slug)
        scenario_slug = 'scenario';
    end
    if startsWith(scenario_slug, 'a0', 'IgnoreCase', true)
        scenario_slug = 'A0_forged';
    elseif startsWith(scenario_slug, 'a1', 'IgnoreCase', true)
        scenario_slug = 'A1_forged';
    elseif startsWith(scenario_slug, 'a2', 'IgnoreCase', true)
        scenario_slug = 'A2_forged';
    elseif contains(scenario_slug, 'baseline')
        scenario_slug = 'baseline_authentic';
    end

    h1 = NaN; h2 = NaN;
    if isfield(rr, 'harmonic_mult_1'), h1 = double(rr.harmonic_mult_1); end
    if isfield(rr, 'harmonic_mult_2'), h2 = double(rr.harmonic_mult_2); end
    if isfinite(h1) && isfinite(h2)
        pair_slug = sprintf('h%d_vs_h%d', round(h1), round(h2));
    else
        pair_slug = sprintf('f%.0f_vs_f%.0f', rr.freq1_hz, rr.freq2_hz);
    end

    out_svg_name = sprintf('D1_%s_%s_enf_consistency', scenario_slug, pair_slug);
    af_save_figure_svg_and_close(f, cfg, out_svg_name);
end


% af_plot_harmonic_consistency_injection_vs_forged - Save a two-panel injection vs forged comparison.
%
%   Plots the injection reference ENF against the forged trace ENF for each
%   of the two harmonic frequencies as side-by-side subplots and saves the figure.
%
%   Inputs:
%     rr_ref    - Consistency result struct for the injection (reference) signal.
%     rr_forged - Consistency result struct for the forged signal.
%     cfg       - AF_CFG struct (for plot export directory resolution).
%
%   Outputs:
%     (none) - Saves D1_<attack>_<pair>_injection_vs_forged_enf_consistency.{png,svg}.
%
function af_plot_harmonic_consistency_injection_vs_forged(rr_ref, rr_forged, cfg)
    if isempty(rr_ref) || isempty(rr_forged)
        return;
    end
    L = min([numel(rr_ref.enf1_fund_hz), numel(rr_ref.enf2_fund_hz), numel(rr_forged.enf1_fund_hz), numel(rr_forged.enf2_fund_hz)]);
    if L <= 1
        return;
    end

    t_ref = rr_ref.time_sec(:);
    t_frg = rr_forged.time_sec(:);
    t_ref = t_ref(1:min(numel(t_ref), L));
    t_frg = t_frg(1:min(numel(t_frg), L));
    L = min([L, numel(t_ref), numel(t_frg)]);
    if L <= 1
        return;
    end
    t = t_frg(1:L);
    e1_ref = rr_ref.enf1_fund_hz(1:L);
    e2_ref = rr_ref.enf2_fund_hz(1:L);
    e1_frg = rr_forged.enf1_fund_hz(1:L);
    e2_frg = rr_forged.enf2_fund_hz(1:L);

    c_h1 = af_calc_pearson_corr(e1_ref, e1_frg);
    c_h2 = af_calc_pearson_corr(e2_ref, e2_frg);

    f = figure('Name', sprintf('%s - Injection vs Forged (D1)', char(rr_forged.scenario_label)), ...
               'Color', 'w', 'Position', [220 120 1100 460]);
    ax1 = subplot(1,2,1, 'Parent', f);
    plot(ax1, t, e1_ref, '-', 'LineWidth', 1.0); hold(ax1, 'on'); grid(ax1, 'on');
    plot(ax1, t, e1_frg, '--', 'LineWidth', 1.0); hold(ax1, 'off');
    xlabel(ax1, 'Time (seconds)');
    ylabel(ax1, 'Implied Fundamental ENF (Hz)');
    title(ax1, sprintf('ENF @ %.1f Hz (implied fundamental) | Corr = %.4f', rr_forged.freq1_hz, c_h1), 'Interpreter', 'none');
    legend(ax1, 'Injection trace', 'Forged trace', 'Location', 'southoutside');

    ax2 = subplot(1,2,2, 'Parent', f);
    plot(ax2, t, e2_ref, '-', 'LineWidth', 1.0); hold(ax2, 'on'); grid(ax2, 'on');
    plot(ax2, t, e2_frg, '--', 'LineWidth', 1.0); hold(ax2, 'off');
    xlabel(ax2, 'Time (seconds)');
    ylabel(ax2, 'Implied Fundamental ENF (Hz)');
    title(ax2, sprintf('ENF @ %.1f Hz (implied fundamental) | Corr = %.4f', rr_forged.freq2_hz, c_h2), 'Interpreter', 'none');
    legend(ax2, 'Injection trace', 'Forged trace', 'Location', 'southoutside');

    try
        sgtitle(f, sprintf('%s | D1 Injection vs Forged ENF overlay', char(rr_forged.scenario_label)), 'Interpreter', 'none');
    catch
        % no-op for older MATLAB
    end

    attack_id = af_attack_id_from_text(rr_forged.scenario_key, false);
    if isempty(attack_id) || strcmpi(attack_id, 'AX')
        attack_id = 'AX';
    end
    h1 = NaN; h2 = NaN;
    if isfield(rr_forged, 'harmonic_mult_1'), h1 = double(rr_forged.harmonic_mult_1); end
    if isfield(rr_forged, 'harmonic_mult_2'), h2 = double(rr_forged.harmonic_mult_2); end
    if isfinite(h1) && isfinite(h2)
        pair_slug = sprintf('h%d_vs_h%d', round(h1), round(h2));
    else
        pair_slug = sprintf('f%.0f_vs_f%.0f', rr_forged.freq1_hz, rr_forged.freq2_hz);
    end
    out_svg_name = sprintf('D1_%s_%s_injection_vs_forged_enf_consistency', lower(af_sanitize_filename(char(attack_id))), pair_slug);
    af_save_figure_svg_and_close(f, cfg, out_svg_name);
end

function T = af_harmonic_consistency_make_table(rows)
    if isempty(rows)
        T = table();
        return;
    end
    n = numel(rows);
    scenario_key = cell(n,1);
    h1_mult = nan(n,1);
    h2_mult = nan(n,1);
    freq1_hz = nan(n,1);
    freq2_hz = nan(n,1);
    corr_norm = nan(n,1);
    for i = 1:n
        r = rows(i);
        scenario_key{i} = char(r.scenario_key);
        if isfield(r, 'harmonic_mult_1'), h1_mult(i) = double(r.harmonic_mult_1); end
        if isfield(r, 'harmonic_mult_2'), h2_mult(i) = double(r.harmonic_mult_2); end
        freq1_hz(i) = double(r.freq1_hz);
        freq2_hz(i) = double(r.freq2_hz);
        corr_norm(i) = double(r.corr_norm);
    end
    T = table(scenario_key, h1_mult, h2_mult, freq1_hz, freq2_hz, corr_norm, ...
              'VariableNames', {'ScenarioKey','H1_Mult','H2_Mult','Freq1_Hz','Freq2_Hz','NormCorr'});
end



function af_harmonic_consistency_print_table(T)
    fprintf('\nINFO[AF][D1_inter_harmonic_consistency_check] Inter-frequency consistency check (normalized ENF overlays, no scatter plot)\n');
    if isempty(T) || height(T) == 0
        fprintf('  No harmonic consistency results.\n');
        return;
    end
    keys = strings(height(T), 1);
    try
        keys = string(T.ScenarioKey);
    catch
        % keep default empty keys
    end

    mask_baseline = strcmpi(keys, "baseline_authentic_trace");
    mask_A2 = startsWith(lower(keys), "a2_2_");
    mask_other = ~(mask_baseline | mask_A2);

    d1_metrics = {'Freq1_Hz','Freq2_Hz','NormCorr'};
    if any(mask_baseline)
        fprintf('  Baseline (Authentic trace) records:\n');
        Tsub = T(mask_baseline, :);
        af_print_table_records_as_columns(Tsub, '', d1_col_headers_from_table(Tsub), d1_metrics);
    end
    if any(mask_A2)
        fprintf('  A2 synthetic embedding variant records:\n');
        Tsub = T(mask_A2, :);
        af_print_table_records_as_columns(Tsub, '', d1_col_headers_from_table(Tsub), d1_metrics);
    end
    if any(mask_other)
        fprintf('  Other scenario records:\n');
        Tsub = T(mask_other, :);
        af_print_table_records_as_columns(Tsub, '', d1_col_headers_from_table(Tsub), d1_metrics);
    end
end


function hdrs = d1_col_headers_from_table(Tsub)
    n = height(Tsub);
    hdrs = cell(1, n);
    for j = 1:n
        h1 = round(double(Tsub.H1_Mult(j)));
        h2 = round(double(Tsub.H2_Mult(j)));
        hdrs{j} = sprintf('H%d_vs_H%d', h1, h2);
    end
end


function af_print_fft_harmonic_comparison_report_combined(x_authentic, ...
        x_case_A0_2, x_case_A1_2, x_case_A2_2, ...
        fs, cfg, runA0_2, runA1_2, runA2_2)
    if nargin < 7, runA0_2 = false; end
    if nargin < 8, runA1_2 = false; end
    if nargin < 9, runA2_2 = false; end
    if nargin < 2, x_case_A0_2 = []; end
    if nargin < 3, x_case_A1_2 = []; end
    if nargin < 4, x_case_A2_2 = []; end

    freq_list = af_get_report_harmonic_freqs(cfg);
    if isempty(freq_list)
        return;
    end

    % ENF (fundamental) row uses peak in +/- enf_bw; all other rows use nearest bin.
    bw_list = zeros(size(freq_list));
    enf_bw = af_get_cfg_numeric(cfg, 'A0_ENF_HALF_BW_HZ', 1.0);
    enf_mask = abs(freq_list - cfg.NOMINAL_FREQ_HZ) < 0.5;
    bw_list(enf_mask) = enf_bw;

    authentic_db = af_fft_harmonic_levels_db(x_authentic, fs, freq_list, bw_list);
    if all(~isfinite(authentic_db))
        fprintf('INFO[AF][Report] Skipping combined FFT harmonic report (authentic trace FFT is empty).\n');
        return;
    end

    includeA0_2 = logical(runA0_2) && ~isempty(x_case_A0_2);
    includeA1_2 = logical(runA1_2) && ~isempty(x_case_A1_2);
    includeA2_2 = logical(runA2_2) && ~isempty(x_case_A2_2);

    a0_2_db = nan(size(freq_list));
    a1_2_db = nan(size(freq_list));
    a2_2_db = nan(size(freq_list));
    if includeA0_2, a0_2_db = af_fft_harmonic_levels_db(x_case_A0_2, fs, freq_list, bw_list); end
    if includeA1_2, a1_2_db = af_fft_harmonic_levels_db(x_case_A1_2, fs, freq_list, bw_list); end
    if includeA2_2, a2_2_db = af_fft_harmonic_levels_db(x_case_A2_2, fs, freq_list, bw_list); end

    fprintf('  Observed FFT harmonic levels (dB of |FFT|; ENF row=peak in +/-%0.1f Hz, others=nearest bin)\n', enf_bw);
    header_parts = cell(1, 4);
    n_header_parts = 1;
    header_parts{n_header_parts} = sprintf('%10s | %22s', 'Freq(Hz)', 'Authentic Trace (dB)');
    if includeA0_2, n_header_parts = n_header_parts + 1; header_parts{n_header_parts} = sprintf('%22s | %15s', 'Forged A0_2 (dB)', 'Atten A0_2 (dB)'); end
    if includeA1_2, n_header_parts = n_header_parts + 1; header_parts{n_header_parts} = sprintf('%22s | %15s', 'Forged A1_2 (dB)', 'Atten A1_2 (dB)'); end
    if includeA2_2, n_header_parts = n_header_parts + 1; header_parts{n_header_parts} = sprintf('%22s | %15s', 'Forged A2_2 (dB)', 'Atten A2_2 (dB)'); end
    header = strjoin(header_parts(1:n_header_parts), ' | ');
    fprintf('  %s\n', header);
    fprintf('  %s\n', repmat('-', 1, numel(header)));

    for i = 1:numel(freq_list)
        row_parts = cell(1, 4);
        n_row_parts = 1;
        row_parts{n_row_parts} = sprintf('%10.1f | %22s', freq_list(i), af_fmt_db_for_report(authentic_db(i)));
        if includeA0_2, n_row_parts = n_row_parts + 1; row_parts{n_row_parts} = sprintf('%22s | %15s', af_fmt_db_for_report(a0_2_db(i)), af_fmt_db_for_report(authentic_db(i) - a0_2_db(i))); end
        if includeA1_2, n_row_parts = n_row_parts + 1; row_parts{n_row_parts} = sprintf('%22s | %15s', af_fmt_db_for_report(a1_2_db(i)), af_fmt_db_for_report(authentic_db(i) - a1_2_db(i))); end
        if includeA2_2, n_row_parts = n_row_parts + 1; row_parts{n_row_parts} = sprintf('%22s | %15s', af_fmt_db_for_report(a2_2_db(i)), af_fmt_db_for_report(authentic_db(i) - a2_2_db(i))); end
        row = strjoin(row_parts(1:n_row_parts), ' | ');
        fprintf('  %s\n', row);
    end
end


% af_fft_harmonic_levels_db - Measure FFT magnitude in dB at a list of frequencies.
%
%   For each frequency in freq_list, returns 20*log10(|FFT|) either at the nearest
%   bin (bw_hz=0) or as the peak within a +/- half-bandwidth window (bw_hz>0).
%
%   Inputs:
%     x        - Signal vector.
%     fs       - Sampling rate in Hz.
%     freq_list - Row vector of frequencies to measure (Hz).
%     bw_hz    - (optional) Per-frequency half-bandwidth row vector; 0 = nearest bin.
%
%   Outputs:
%     db_vals - Row vector of dB values corresponding to freq_list.
%
function db_vals = af_fft_harmonic_levels_db(x, fs, freq_list, bw_hz)
    % bw_hz: optional per-frequency half-bandwidth.
    %        0 = nearest-bin (default). >0 = max(|FFT|) in [f-bw, f+bw].
    if nargin < 4 || isempty(bw_hz)
        bw_hz = zeros(size(freq_list));
    elseif isscalar(bw_hz)
        bw_hz = repmat(double(bw_hz), size(freq_list));
    end
    db_vals = nan(size(freq_list));
    [f_fft, m_fft] = af_fft_mag_abs_no_dc(x, fs);
    if isempty(f_fft) || isempty(m_fft)
        return;
    end
    for i = 1:numel(freq_list)
        bw = max(0, bw_hz(i));
        if bw > 0
            idx = (f_fft >= (freq_list(i) - bw)) & (f_fft <= (freq_list(i) + bw));
            if any(idx)
                db_vals(i) = 20 * log10(max(m_fft(idx)) + eps);
                continue;
            end
        end
        [~, idx0] = min(abs(f_fft - freq_list(i)));
        db_vals(i) = 20 * log10(m_fft(idx0) + eps);
    end
end


function s = af_fmt_db_for_report(v)
    if isfinite(v)
        s = sprintf('%.2f', double(v));
    else
        s = 'N/A';
    end
end


% af_print_a0_bandstop_filter_design_spec - Print the bandstop filter design specification.
%
%   For each harmonic multiplier, designs the equiripple bandstop filter, computes
%   the single-pass and cascade-equivalent attenuation at the stopband centre and
%   worst-case stopband edge, and prints the results in a formatted table.
%
%   Inputs:
%     fs  - Sampling rate in Hz.
%     cfg - AF_CFG struct with A0_ENF_HALF_BW_HZ, A0_TRANSITION_BW_HZ,
%           A0_FIR_ORDER, A0_FIR_WEIGHTS, A0_REPORT_FILTER_RESP_NFFT,
%           ATTACK_HARMONIC_MULTS, NOMINAL_FREQ_HZ.
%
%   Outputs:
%     (none) - Prints filter spec table to console.
%
function af_print_a0_bandstop_filter_design_spec(fs, cfg)
    mults = cfg.ATTACK_HARMONIC_MULTS(:).';
    mults = unique(mults(isfinite(mults) & mults > 0), 'stable');
    if isempty(mults)
        return;
    end

    use_vis = af_get_cfg_logical(cfg, 'A_USE_VISIBLE_FILTERS', false);
    if use_vis
        half_bw  = af_get_cfg_numeric(cfg, 'A_VIS_ENF_HALF_BW_HZ',    5.0);
        trans_bw = af_get_cfg_numeric(cfg, 'A_VIS_TRANSITION_BW_HZ',   1.0);
        fir_ord  = round(af_get_cfg_numeric(cfg, 'A_VIS_FIR_ORDER',    2700));
        W = [1 1 1];
        if isfield(cfg, 'A_VIS_FIR_WEIGHTS') && isnumeric(cfg.A_VIS_FIR_WEIGHTS) && numel(cfg.A_VIS_FIR_WEIGHTS) == 3
            W = double(cfg.A_VIS_FIR_WEIGHTS(:).');
        end
        profile_tag = ' [VIS wide-band profile]';
    else
        half_bw  = cfg.A0_ENF_HALF_BW_HZ;
        trans_bw = cfg.A0_TRANSITION_BW_HZ;
        fir_ord  = round(cfg.A0_FIR_ORDER);
        W = [1 1 1];
        if isfield(cfg, 'A0_FIR_WEIGHTS') && isnumeric(cfg.A0_FIR_WEIGHTS) && numel(cfg.A0_FIR_WEIGHTS) == 3
            W = double(cfg.A0_FIR_WEIGHTS(:).');
        end
        profile_tag = '';
    end

    nresp = 262144;
    if isfield(cfg, 'A0_REPORT_FILTER_RESP_NFFT') && isfinite(cfg.A0_REPORT_FILTER_RESP_NFFT) && cfg.A0_REPORT_FILTER_RESP_NFFT > 1024
        nresp = round(cfg.A0_REPORT_FILTER_RESP_NFFT);
    end

    fprintf('  Analysis Fs (Hz): %.2f%s\n', fs, profile_tag);
    fprintf('  Nominal ENF (Hz): %.2f\n', cfg.NOMINAL_FREQ_HZ);
    fprintf('  Harmonic multipliers: ');
    fprintf('%g ', mults);
    fprintf('\n');
    fprintf('  ENF half-bandwidth (Hz): %.2f\n', half_bw);
    fprintf('  Transition bandwidth (Hz): %.2f\n', trans_bw);
    fprintf('  FIR order: %d\n', fir_ord);
    fprintf('  FIR weights [left stop, pass, right stop]: [%.3g %.3g %.3g]\n', W(1), W(2), W(3));
    fprintf('  Stopband per harmonic: [f0-%.2f, f0+%.2f] Hz\n', half_bw, half_bw);
    fprintf('  Note: single-pass compensated FIR application (H(f), not H^2(f)).\n\n');

    fprintf('  %5s | %8s | %19s | %19s | %14s | %14s | %14s | %14s\n', ...
            'Harm', 'f0(Hz)', 'Left Transition', 'Notch Band', ...
            'Center 1x(dB)', 'Center 2x(dB)', 'Worst 1x(dB)', 'Worst 2x(dB)');
    fprintf('  %s\n', repmat('-', 1, 133));

    for m = mults
        f0 = m * cfg.NOMINAL_FREQ_HZ;
        nyq = fs / 2;
        f_pass_1 = f0 - half_bw;
        f_pass_2 = f0 + half_bw;
        f_stop_1 = f_pass_1 - trans_bw;
        f_stop_2 = f_pass_2 + trans_bw;
        ok_edges = (f_stop_1 > 0) && (f_stop_2 < nyq) && (f_stop_1 < f_pass_1) && (f_pass_2 < f_stop_2);
        if ok_edges
            left_tr = sprintf('%.2f-%.2f', f_stop_1, f_pass_1);
            notch_tr = sprintf('%.2f-%.2f', f_pass_1, f_pass_2);
        else
            left_tr = 'N/A';
            notch_tr = 'N/A';
        end

        [~, b_bs, is_valid] = af_design_paper_equiripple_filters(fs, f0, cfg);
        if ~is_valid || isempty(b_bs)
            fprintf('  %5.0f | %8.1f | %19s | %19s | %14s | %14s | %14s | %14s\n', ...
                    m, f0, left_tr, notch_tr, 'N/A', 'N/A', 'N/A', 'N/A');
            continue;
        end

        [f_grid, mag_lin] = af_fir_mag_response_fft(b_bs, fs, nresp);
        if isempty(f_grid) || isempty(mag_lin)
            fprintf('  %5.0f | %8.1f | %19s | %19s | %14s | %14s | %14s | %14s\n', ...
                    m, f0, left_tr, notch_tr, 'N/A', 'N/A', 'N/A', 'N/A');
            continue;
        end

        f_sb1 = f0 - half_bw;
        f_sb2 = f0 + half_bw;
        idx_sb = (f_grid >= f_sb1) & (f_grid <= f_sb2);
        if ~any(idx_sb)
            fprintf('  %5.0f | %8.1f | %19s | %19s | %14s | %14s | %14s | %14s\n', ...
                    m, f0, left_tr, notch_tr, 'N/A', 'N/A', 'N/A', 'N/A');
            continue;
        end

        [~, idx_c] = min(abs(f_grid - f0));
        Hc = mag_lin(idx_c);
        Hsb = mag_lin(idx_sb);

        % Positive attenuation numbers in dB.
        att_center_1x = -20 * log10(Hc + eps);
        att_center_2x = -20 * log10((Hc.^2) + eps);  % effective for forward-backward

        % Worst attenuation = least attenuation across stopband (largest gain in stopband).
        Hsb_max = max(Hsb);
        att_worst_1x = -20 * log10(Hsb_max + eps);
        att_worst_2x = -20 * log10((Hsb_max.^2) + eps);

        fprintf('  %5.0f | %8.1f | %19s | %19s | %14.2f | %14.2f | %14.2f | %14.2f\n', ...
                m, f0, left_tr, notch_tr, att_center_1x, att_center_2x, att_worst_1x, att_worst_2x);
    end
end


% af_get_report_harmonic_freqs - Return the list of harmonic frequencies for FFT reporting.
%
%   Reads from FFT_HARMONIC_TEXT_FREQS_HZ, falls back to ATTACK_HARMONIC_MULTS x
%   NOMINAL_FREQ_HZ, and further falls back to the default 60 Hz series if neither
%   is configured. Filters to positive finite values not exceeding 500 Hz.
%
%   Inputs:
%     cfg - AF_CFG struct.
%
%   Outputs:
%     freq_list - Sorted unique row vector of harmonic frequencies (Hz).
%
function freq_list = af_get_report_harmonic_freqs(cfg)
    if isfield(cfg, 'FFT_HARMONIC_TEXT_FREQS_HZ') && ~isempty(cfg.FFT_HARMONIC_TEXT_FREQS_HZ)
        freq_list = cfg.FFT_HARMONIC_TEXT_FREQS_HZ(:).';
    elseif isfield(cfg, 'ATTACK_HARMONIC_MULTS') && isfield(cfg, 'NOMINAL_FREQ_HZ')
        freq_list = cfg.ATTACK_HARMONIC_MULTS(:).' * cfg.NOMINAL_FREQ_HZ;
    else
        freq_list = [60 120 180 240 300 360 420];
    end
    freq_list = unique(freq_list(isfinite(freq_list) & freq_list > 0), 'stable');
    freq_list = freq_list(freq_list <= 500);
end


% af_fir_mag_response_fft - Compute the one-sided FFT magnitude response of an FIR filter.
%
%   Inputs:
%     b        - FIR filter coefficient vector.
%     fs       - Sampling rate in Hz (used to build the frequency axis).
%     nfft_resp - FFT size for the response computation.
%
%   Outputs:
%     f_grid  - One-sided frequency axis vector (Hz), from 0 to fs/2.
%     mag_lin - Linear magnitude vector corresponding to f_grid.
%
function [f_grid, mag_lin] = af_fir_mag_response_fft(b, fs, nfft_resp)
    b = b(:);
    if isempty(b) || nfft_resp < 8
        f_grid = [];
        mag_lin = [];
        return;
    end
    H = fft(b, nfft_resp);
    H = H(1:floor(nfft_resp/2)+1);
    mag_lin = abs(H);
    f_grid = (0:floor(nfft_resp/2))' * (fs / nfft_resp);
end


% af_load_mono_trace_resampled - Load a WAV file as mono, resample, and optionally normalize.
%
%   Reads the WAV file at wav_path, converts to mono, removes DC, replaces
%   non-finite samples with zero, resamples to target_fs, and optionally applies
%   unit-peak normalisation.
%
%   Inputs:
%     wav_path              - Full path to the input WAV file.
%     target_fs             - Target sampling rate in Hz.
%     do_unit_peak_normalize - (optional) Logical; if true, divide by peak abs value
%                             (default: true).
%
%   Outputs:
%     x_out  - Resampled mono signal vector at target_fs Hz.
%     fs_out - Output sampling rate (equals target_fs).
%
function [x_out, fs_out] = af_load_mono_trace_resampled(wav_path, target_fs, do_unit_peak_normalize)
    if nargin < 3 || isempty(do_unit_peak_normalize)
        do_unit_peak_normalize = true;
    end
    if ~exist(wav_path, 'file')
        error('ERROR: AF input file not found: %s', wav_path);
    end

    [x, fs_in] = audioread(wav_path);
    if isempty(x)
        x_out = [];
        fs_out = target_fs;
        return;
    end

    % Mono conversion for consistent ENF-band processing.
    if size(x, 2) > 1
        x = mean(x, 2);
    else
        x = x(:, 1);
    end

    x = x(:);
    x = x - mean(x, 'omitnan');
    x(~isfinite(x)) = 0;

    if fs_in == target_fs
        x_out = x;
    else
        x_out = resample(x, target_fs, fs_in);
    end

    x_out = x_out(:);
    if do_unit_peak_normalize && any(isfinite(x_out))
        mx = max(abs(x_out));
        if isfinite(mx) && mx > 0
            x_out = x_out / mx;
        end
    end

    fs_out = target_fs;
end


% af_apply_enf_bandstop_removal - Apply A0 equiripple bandstop FIR per harmonic.
%
%   Applies a single-pass compensated equiripple bandstop FIR filter to each
%   configured harmonic multiplier to realize the designed H(f) attenuation
%   without the H(f)^2 squaring that forward-backward filtering would produce.
%
%   Inputs:
%     x   - Input signal vector.
%     fs  - Sampling rate in Hz.
%     cfg - AF_CFG struct with ATTACK_HARMONIC_MULTS, NOMINAL_FREQ_HZ, and
%           all A0 filter design parameters.
%
%   Outputs:
%     x_removed - Signal with ENF band attenuated at each configured harmonic.
%
function x_removed = af_apply_enf_bandstop_removal(x, fs, cfg)
    % A0 removal path: single-pass compensated FIR per harmonic to realize
    % designed H(f) attenuation (not H(f)^2 from forward-backward filtering).
    x_removed = x(:);
    mults = cfg.ATTACK_HARMONIC_MULTS(:).';
    for m = mults
        center_hz = m * cfg.NOMINAL_FREQ_HZ;
        [~, b_bs, is_valid] = af_design_paper_equiripple_filters(fs, center_hz, cfg);
        if ~is_valid
            warning('AF:skipBandstop', 'Skipping invalid harmonic band for %g Hz (outside realizable range).', center_hz);
            continue;
        end
        x_removed = af_apply_fir_single_pass_compensated(x_removed, b_bs);
    end
    x_removed = x_removed(:);
end


% af_apply_fir_single_pass_compensated - Apply a linear-phase FIR with group-delay compensation.
%
%   Applies the filter in the forward direction only, then advances the output
%   by order/2 samples to compensate the constant group delay of a symmetric
%   even-order FIR. Preserves H(f) rather than H(f)^2 as filtfilt would.
%
%   Inputs:
%     x - Input signal vector.
%     b - FIR filter coefficient vector (even order required).
%
%   Outputs:
%     y - Filtered signal with group delay removed (same length as x).
%
function y = af_apply_fir_single_pass_compensated(x, b)
    % Single-pass FIR with group-delay compensation:
    % preserves the designed H(f) magnitude response (no squaring as in filtfilt).
    x = x(:);
    b = b(:);
    if isempty(x) || isempty(b)
        y = x;
        return;
    end

    order = numel(b) - 1;
    if mod(order, 2) ~= 0
        warning('AF:OddOrderFallback', ...
                ['FIR order %d is odd; non-integer group delay. ' ...
                 'Falling back to zero-phase application.'], order);
        y = af_apply_fir_zero_phase(x, b);
        return;
    end

    delay = order / 2;
    if numel(x) <= delay
        y = zeros(size(x));
        return;
    end

    y_delayed = filter(b, 1, x);
    y = [y_delayed(delay+1:end); zeros(delay, 1)];
    y = y(:);
end

% af_build_a1_from_a0_variant - Produce the A1 noise-fill signal on top of an A0 base.
%
%   Generates the narrow-band noise fill matched to the authentic spectrum and
%   adds it to x_a0_base to produce the A1 forged signal.
%
%   Inputs:
%     x_authentic - Original authentic signal vector.
%     x_a0_base   - A0 ENF-removed signal to use as the A1 base.
%     fs          - Sampling rate in Hz.
%     cfg_noise   - AF_CFG struct with A1 noise generation parameters.
%
%   Outputs:
%     x_case_A1  - A1 forged signal (A0 base + noise fill).
%     x_noise_fill - The bandpass noise fill component alone.
%
function [x_case_A1, x_noise_fill] = af_build_a1_from_a0_variant(x_authentic, x_a0_base, fs, cfg_noise)
    x_authentic = x_authentic(:);
    x_a0_base = x_a0_base(:);

    n_samples = min(numel(x_authentic), numel(x_a0_base));
    x_authentic = x_authentic(1:n_samples);
    x_a0_base = x_a0_base(1:n_samples);

    x_noise_fill = af_generate_a1_bandpassed_noise_fill(x_authentic, n_samples, fs, cfg_noise);
    x_case_A1 = x_a0_base + x_noise_fill;
    x_case_A1 = x_case_A1(:);
end


% removed deprecated function af_build_a1_peak_matched_embedding (unused)



% af_generate_a1_bandpassed_noise_fill - Generate spectrally matched narrow-band noise fill.
%
%   For each harmonic multiplier, generates white noise, shapes it with the
%   equiripple bandpass FIR (single-pass compensated, consistent with A0),
%   and scales it so its in-band level matches the median of the authentic
%   signal's neighbouring spectrum.
%
%   Inputs:
%     x_ref    - Authentic signal vector used as the spectral reference.
%     n_samples - Number of samples to generate.
%     fs       - Sampling rate in Hz.
%     cfg      - AF_CFG struct with A1 and A0 noise/filter parameters.
%
%   Outputs:
%     x_noise_fill - Bandpassed noise fill vector (n_samples x 1).
%
function x_noise_fill = af_generate_a1_bandpassed_noise_fill(x_ref, n_samples, fs, cfg)
    % A1/A2 attack-generation convention:
    % bandpass shaping is single-pass compensated FIR (H(f), not H(f)^2).
    % This keeps insertion/removal shaping consistent with A0.
    x_ref = x_ref(:);
    x_noise_fill = zeros(n_samples, 1);
    mults = cfg.ATTACK_HARMONIC_MULTS(:).';
    seed0 = round(af_get_cfg_numeric(cfg, 'A1_RANDOM_SEED', 42));
    stat_mode = 'median';
    if isfield(cfg, 'A1_NOISE_MATCH_STAT') && ~isempty(cfg.A1_NOISE_MATCH_STAT)
        stat_mode = lower(char(cfg.A1_NOISE_MATCH_STAT));
    end

    for ii = 1:numel(mults)
        m = mults(ii);
        center_hz = m * cfg.NOMINAL_FREQ_HZ;
        [b_bp, ~, is_valid] = af_design_paper_equiripple_filters(fs, center_hz, cfg);
        if ~is_valid || isempty(b_bp)
            continue;
        end

        % Reproducible but different noise for each harmonic.
        try
            rng(seed0 + 100 + ii);
        catch
        end
        n = randn(n_samples, 1);
        n_bp = af_apply_fir_single_pass_compensated(n, b_bp);

        [f_pass_1, f_pass_2, f_stop_1, f_stop_2, ok_edges] = af_get_a1_band_edges(center_hz, fs, cfg);
        if ~ok_edges
            x_noise_fill = x_noise_fill + n_bp;
            continue;
        end

        if af_get_cfg_logical(cfg, 'A_USE_VISIBLE_FILTERS', false)
            default_neigh_bw = 2 * cfg.A_VIS_ENF_HALF_BW_HZ;
        else
            default_neigh_bw = 2 * cfg.A0_ENF_HALF_BW_HZ;
        end
        neigh_bw = af_get_cfg_numeric(cfg, 'A1_NOISE_NEIGHBOR_BW_HZ', default_neigh_bw);
        neigh_bw = max(eps, neigh_bw);

        bands_neigh = [max(0, f_stop_1 - neigh_bw), f_stop_1; ...
                       f_stop_2, min(fs/2, f_stop_2 + neigh_bw)];
        bands_neigh = af_prune_invalid_bands(bands_neigh);
        band_enf = [f_pass_1, f_pass_2];

        target_level = af_fft_band_stat_abs(x_ref, fs, bands_neigh, stat_mode);
        src_level = af_fft_band_stat_abs(n_bp, fs, band_enf, stat_mode);

        if isfinite(target_level) && isfinite(src_level) && src_level > 0
            n_bp = n_bp * (target_level / src_level);
        end
        x_noise_fill = x_noise_fill + n_bp(:);
    end

    x_noise_fill = x_noise_fill(:);
end


% af_a1_match_syn_peak_magnitude - Find the scaling factor M so A2 peak matches the authentic.
%
%   Uses a bounded binary search to find M such that the FFT peak of
%   (x_base + M * x_syn_unit) within the configured bandwidth around the nominal
%   ENF frequency equals the peak of the authentic trace. Returns M=0 if the base
%   peak already meets or exceeds the target.
%
%   Inputs:
%     x_authentic - Authentic signal vector (used to measure the target peak).
%     x_base      - A1 base signal (x_case_A1_2).
%     x_syn_unit  - Unit synthetic ENF component vector.
%     fs          - Sampling rate in Hz.
%     cfg         - AF_CFG struct with A2_PEAK_MATCH_FREQ_HZ, A2_PEAK_MATCH_BW_HZ.
%
%   Outputs:
%     M_opt - Optimal scaling factor M (non-negative scalar).
%     info  - Struct with fields: peak_match_freq_hz, peak_match_bw_hz,
%             target_peak_abs, base_peak_abs, final_peak_abs,
%             binsearch_iters, tol_abs.
%
function [M_opt, info] = af_a1_match_syn_peak_magnitude(x_authentic, x_base, x_syn_unit, fs, cfg)
    peak_freq_hz = af_get_cfg_numeric(cfg, 'A2_PEAK_MATCH_FREQ_HZ', af_get_cfg_numeric(cfg, 'A1_PEAK_MATCH_FREQ_HZ', cfg.NOMINAL_FREQ_HZ));
    peak_bw_hz = af_get_cfg_numeric(cfg, 'A2_PEAK_MATCH_BW_HZ', af_get_cfg_numeric(cfg, 'A1_PEAK_MATCH_BW_HZ', cfg.A0_ENF_HALF_BW_HZ));
    peak_bw_hz = max(0, peak_bw_hz);
    n_iter = 20;
    n_iter = max(1, n_iter);
    tol_frac = 1e-3;
    tol_frac = max(0, tol_frac);

    target_peak = af_fft_local_peak_abs(x_authentic, fs, peak_freq_hz, peak_bw_hz);
    base_peak = af_fft_local_peak_abs(x_base, fs, peak_freq_hz, peak_bw_hz);

    if ~isfinite(target_peak) || target_peak <= 0
        target_peak = max(target_peak, eps);
    end
    if ~isfinite(base_peak) || base_peak < 0
        base_peak = 0;
    end

    tol_abs = tol_frac * max(target_peak, eps);
    M_lo = 0;
    M_hi = 1;
    p_lo = base_peak;
    p_hi = af_fft_local_peak_abs(x_base + M_hi*x_syn_unit, fs, peak_freq_hz, peak_bw_hz);

    % Expand upper bound until we bracket the target or hit a practical limit.
    expand_count = 0;
    while isfinite(p_hi) && (p_hi < target_peak) && expand_count < 24
        M_hi = 2 * M_hi;
        p_hi = af_fft_local_peak_abs(x_base + M_hi*x_syn_unit, fs, peak_freq_hz, peak_bw_hz);
        expand_count = expand_count + 1;
    end

    if ~isfinite(p_hi)
        p_hi = p_lo;
    end

    % If base already exceeds target, keep M = 0 (noise fill alone is enough/high).
    if p_lo >= target_peak
        M_opt = 0;
        final_peak = p_lo;
    else
        M_best = M_hi;
        p_best = p_hi;
        if ~isfinite(p_best)
            M_best = 0;
            p_best = p_lo;
        end

        for k = 1:n_iter
            M_mid = 0.5 * (M_lo + M_hi);
            p_mid = af_fft_local_peak_abs(x_base + M_mid*x_syn_unit, fs, peak_freq_hz, peak_bw_hz);
            if ~isfinite(p_mid)
                p_mid = p_lo;
            end

            if abs(p_mid - target_peak) < abs(p_best - target_peak)
                M_best = M_mid;
                p_best = p_mid;
            end
            if abs(p_mid - target_peak) <= tol_abs
                M_best = M_mid;
                p_best = p_mid;
                break;
            end

            if p_mid < target_peak
                M_lo = M_mid;
                p_lo = p_mid;
            else
                M_hi = M_mid;
            end
        end

        M_opt = max(0, M_best);
        final_peak = p_best;
    end

    info = struct();
    info.peak_match_freq_hz = peak_freq_hz;
    info.peak_match_bw_hz = peak_bw_hz;
    info.target_peak_abs = target_peak;
    info.base_peak_abs = base_peak;
    info.final_peak_abs = final_peak;
    info.binsearch_iters = n_iter;
    info.tol_abs = tol_abs;
end


% af_fft_local_peak_abs - Return the peak one-sided FFT magnitude near a target frequency.
%
%   Finds the maximum |FFT| within center_hz +/- half_bw_hz, or at the nearest
%   bin if half_bw_hz is zero or absent.
%
%   Inputs:
%     x           - Signal vector.
%     fs          - Sampling rate in Hz.
%     center_hz   - Centre frequency to search around (Hz).
%     half_bw_hz  - (optional) Half-bandwidth of the search window (Hz); 0 = nearest bin.
%
%   Outputs:
%     peak_abs - Maximum linear |FFT| in the specified window, or NaN if unavailable.
%
function peak_abs = af_fft_local_peak_abs(x, fs, center_hz, half_bw_hz)
    [f_fft, mag_abs] = af_fft_mag_abs_no_dc(x, fs);
    if isempty(f_fft) || isempty(mag_abs)
        peak_abs = NaN;
        return;
    end

    if nargin < 4 || ~isfinite(half_bw_hz) || half_bw_hz < 0
        half_bw_hz = 0;
    end

    if half_bw_hz > 0
        idx = (f_fft >= (center_hz - half_bw_hz)) & (f_fft <= (center_hz + half_bw_hz));
        if any(idx)
            peak_abs = max(mag_abs(idx));
            return;
        end
    end

    [~, idx0] = min(abs(f_fft - center_hz));
    peak_abs = mag_abs(idx0);
end


% af_fft_band_stat_abs - Compute a statistic of |FFT| magnitude over one or more bands.
%
%   Collects all FFT magnitude values within the union of the specified frequency
%   bands and returns their median (default) or mean.
%
%   Inputs:
%     x        - Signal vector.
%     fs       - Sampling rate in Hz.
%     bands_hz - Nx2 matrix where each row is [f_lo, f_hi] in Hz.
%     stat_mode - (optional) 'median' (default) or 'mean'.
%
%   Outputs:
%     val - Scalar statistic of |FFT| within the union of the bands, or NaN.
%
function val = af_fft_band_stat_abs(x, fs, bands_hz, stat_mode)
    [f_fft, mag_abs] = af_fft_mag_abs_no_dc(x, fs);
    if isempty(f_fft) || isempty(mag_abs) || isempty(bands_hz)
        val = NaN;
        return;
    end

    if nargin < 4 || isempty(stat_mode)
        stat_mode = 'median';
    end

    mask = false(size(f_fft));
    for k = 1:size(bands_hz, 1)
        b1 = bands_hz(k,1);
        b2 = bands_hz(k,2);
        if ~isfinite(b1) || ~isfinite(b2) || b2 <= b1
            continue;
        end
        mask = mask | ((f_fft >= b1) & (f_fft <= b2));
    end

    vals = mag_abs(mask);
    vals = vals(isfinite(vals));
    if isempty(vals)
        val = NaN;
        return;
    end
    val = af_summarize_values(vals, stat_mode);
end


% af_prune_invalid_bands - Remove band rows where the lower edge meets or exceeds the upper.
%
%   Inputs:
%     bands_in - Nx2 matrix of [f_lo, f_hi] band edges.
%
%   Outputs:
%     bands_out - Subset of rows where f_lo < f_hi and both are finite.
%
function bands_out = af_prune_invalid_bands(bands_in)
    bands_out = [];
    if isempty(bands_in)
        return;
    end
    for k = 1:size(bands_in,1)
        b1 = bands_in(k,1);
        b2 = bands_in(k,2);
        if isfinite(b1) && isfinite(b2) && (b2 > b1)
            bands_out = [bands_out; b1 b2]; %#ok<AGROW>
        end
    end
end


% af_get_a1_band_edges - Compute the A1 bandpass filter passband and stopband edges.
%
%   Returns the four edge frequencies for the equiripple bandpass filter centred
%   at center_hz, derived from the A0 (or visible variant) half-bandwidth and
%   transition bandwidth parameters.
%
%   Inputs:
%     center_hz - Centre frequency of the ENF band (Hz).
%     fs        - Sampling rate in Hz.
%     cfg       - AF_CFG struct with A0_ENF_HALF_BW_HZ and A0_TRANSITION_BW_HZ
%                 (or A_VIS_* equivalents if A_USE_VISIBLE_FILTERS is true).
%
%   Outputs:
%     f_pass_1 - Lower passband edge (Hz).
%     f_pass_2 - Upper passband edge (Hz).
%     f_stop_1 - Lower stopband edge (Hz).
%     f_stop_2 - Upper stopband edge (Hz).
%     ok       - Logical; true if all edges are within (0, fs/2) and ordered correctly.
%
function [f_pass_1, f_pass_2, f_stop_1, f_stop_2, ok] = af_get_a1_band_edges(center_hz, fs, cfg)
    nyq = fs / 2;
    if af_get_cfg_logical(cfg, 'A_USE_VISIBLE_FILTERS', false)
        half_bw = cfg.A_VIS_ENF_HALF_BW_HZ;
        trans_bw = cfg.A_VIS_TRANSITION_BW_HZ;
    else
        half_bw = cfg.A0_ENF_HALF_BW_HZ;
        trans_bw = cfg.A0_TRANSITION_BW_HZ;
    end
    f_pass_1 = center_hz - half_bw;
    f_pass_2 = center_hz + half_bw;
    f_stop_1 = f_pass_1 - trans_bw;
    f_stop_2 = f_pass_2 + trans_bw;
    ok = isfinite(f_pass_1) && isfinite(f_pass_2) && isfinite(f_stop_1) && isfinite(f_stop_2) && ...
         (f_stop_1 > 0) && (f_stop_2 < nyq) && (f_stop_1 < f_pass_1) && (f_pass_2 < f_stop_2);
end


% af_generate_syn_enf_component - Generate the A2 synthetic ENF injection component.
%
%   Builds a smooth frequency-modulated trajectory around NOMINAL_FREQ_HZ using
%   sinusoidal modulation terms and a low-amplitude smoothed-noise drift, then
%   synthesizes each configured harmonic and projects it through the paper-spec
%   equiripple bandpass FIR (single-pass compensated, consistent with A0/A1).
%
%   Inputs:
%     n_samples - Number of output samples to generate.
%     fs        - Sampling rate in Hz.
%     cfg       - AF_CFG struct with A2_SYN_MOD_FREQS_HZ, A2_SYN_MOD_AMPS_HZ,
%                 A2_SYN_ENF_DEV_HZ, A2_RANDOM_SEED, ATTACK_HARMONIC_MULTS,
%                 NOMINAL_FREQ_HZ, and all A0 filter design parameters.
%
%   Outputs:
%     x_syn - Synthetic ENF component vector (n_samples x 1).
%
function x_syn = af_generate_syn_enf_component(n_samples, fs, cfg)
    % Generate a smooth syn ENF trajectory (fundamental) and synthesize
    % requested harmonics, then project each harmonic through the same
    % paper-spec bandpass to keep the injected component confined.
    % Attack-generation convention: single-pass compensated FIR shaping
    % (H(f), not H(f)^2) for consistency with A0/A1.
    seed = round(af_get_cfg_numeric(cfg, 'A2_RANDOM_SEED', af_get_cfg_numeric(cfg, 'A1_RANDOM_SEED', 42)));
    rng(seed);

    t = (0:n_samples-1).' / fs;
    delta_f = zeros(n_samples, 1);

    mod_freqs = cfg.A2_SYN_MOD_FREQS_HZ(:);
    mod_amps  = cfg.A2_SYN_MOD_AMPS_HZ(:);
    n_terms = min(numel(mod_freqs), numel(mod_amps));

    for k = 1:n_terms
        phi0 = 2*pi*rand();
        delta_f = delta_f + mod_amps(k) * sin(2*pi*mod_freqs(k)*t + phi0);
    end

    % Add a low-amplitude smoothed noise term for non-periodic drift.
    noise = randn(n_samples, 1);
    smooth_len = max(5, round(2 * fs));  % ~2 s moving-average smoothing
    smooth_kernel = ones(smooth_len, 1) / smooth_len;
    drift = conv(noise, smooth_kernel, 'same');
    drift = drift / max(std(drift, 1, 'omitnan'), eps);
    delta_f = delta_f + 0.02 * drift;

    % Constrain max deviation to the configured bound.
    max_dev = max(abs(delta_f));
    if isfinite(max_dev) && max_dev > 0
        max_dev_hz = af_get_cfg_numeric(cfg, 'A2_SYN_ENF_DEV_HZ', af_get_cfg_numeric(cfg, 'A1_SYN_ENF_DEV_HZ', 0.15));
        delta_f = (max_dev_hz / max_dev) * delta_f;
    end

    f_syn_fund = cfg.NOMINAL_FREQ_HZ + delta_f;
    x_syn = zeros(n_samples, 1);

    mults = cfg.ATTACK_HARMONIC_MULTS(:).';
    for m = mults
        inst_freq = m * f_syn_fund;
        phi = 2*pi*cumsum(inst_freq) / fs + 2*pi*rand();
        y = sin(phi);

        center_hz = m * cfg.NOMINAL_FREQ_HZ;
        [b_bp, ~, is_valid] = af_design_paper_equiripple_filters(fs, center_hz, cfg);
        if is_valid
            y = af_apply_fir_single_pass_compensated(y, b_bp);
        end

        x_syn = x_syn + y(:);
    end

    x_syn = detrend(x_syn - mean(x_syn, 'omitnan'));
end


% af_design_paper_equiripple_filters - Design the paper-spec equiripple bandpass and bandstop FIR pair.
%
%   Constructs the six-point frequency grid centred at center_hz using A0_ENF_HALF_BW_HZ
%   and A0_TRANSITION_BW_HZ, then calls firpm (preferred), remez (fallback), or a
%   windowed-sinc approximation (last resort). Returns both a bandpass and a complementary
%   bandstop filter. Delegates to af_design_visible_variant_filters when A_USE_VISIBLE_FILTERS
%   is set in cfg.
%
%   Inputs:
%     fs        - Sampling rate in Hz.
%     center_hz - Centre frequency of the ENF band (Hz).
%     cfg       - AF_CFG struct with A0_ENF_HALF_BW_HZ, A0_TRANSITION_BW_HZ,
%                 A0_FIR_ORDER, A0_FIR_WEIGHTS, and A_USE_VISIBLE_FILTERS.
%
%   Outputs:
%     b_bp     - Bandpass FIR coefficient vector.
%     b_bs     - Bandstop FIR coefficient vector.
%     is_valid - Logical; false if the requested band falls outside the realizable range.
%
function [b_bp, b_bs, is_valid] = af_design_paper_equiripple_filters(fs, center_hz, cfg)
    % Equiripple filter layout for ENF in-band and transition regions:
    %   ENF in-band: center +/- A0_ENF_HALF_BW_HZ
    %   Transition width: A0_TRANSITION_BW_HZ on each side
    % Example (60 Hz): stop/pass edges at 51, 59, 61, 69 Hz
    if af_get_cfg_logical(cfg, 'A_USE_VISIBLE_FILTERS', false)
        [b_bp, b_bs, is_valid] = af_design_visible_variant_filters(fs, center_hz, cfg);
        return;
    end

    nyq = fs / 2;
    persistent af_warned_remez_fallback af_warned_windowed_fallback

    f_pass_1 = center_hz - cfg.A0_ENF_HALF_BW_HZ;
    f_pass_2 = center_hz + cfg.A0_ENF_HALF_BW_HZ;
    f_stop_1 = f_pass_1 - cfg.A0_TRANSITION_BW_HZ;
    f_stop_2 = f_pass_2 + cfg.A0_TRANSITION_BW_HZ;

    is_valid = true;
    if f_stop_1 <= 0 || f_stop_2 >= nyq || f_pass_1 <= 0 || f_pass_2 >= nyq || ~(f_stop_1 < f_pass_1 && f_pass_2 < f_stop_2)
        is_valid = false;
        b_bp = [];
        b_bs = [];
        return;
    end

    N = cfg.A0_FIR_ORDER;
    W = cfg.A0_FIR_WEIGHTS;
    if numel(W) ~= 3
        W = [1 1 1];
    end

    f_edges = [0, f_stop_1, f_pass_1, f_pass_2, f_stop_2, nyq] / nyq;
    a_bp = [0 0 1 1 0 0];
    a_bs = [1 1 0 0 1 1];

    % Preferred: paper-faithful equiripple design via firpm/remez.
    firpm_exist_code = exist('firpm', 'file');
    if ismember(firpm_exist_code, [2 3 5 6])
        b_bp = firpm(N, f_edges, a_bp, W);
        b_bs = firpm(N, f_edges, a_bs, W);
        return;
    end
    remez_exist_code = exist('remez', 'file');
    if ismember(remez_exist_code, [2 3 5 6])
        b_bp = remez(N, f_edges, a_bp, W);
        b_bs = remez(N, f_edges, a_bs, W);
        if isempty(af_warned_remez_fallback) || ~af_warned_remez_fallback
            warning('AF:firpmUnavailableUsingRemez', ...
                    'firpm not found. Using remez() as equiripple fallback for anti-forensics filter design.');
            af_warned_remez_fallback = true;
        end
        return;
    end

    % Last-resort fallback (no equiripple available): windowed-sinc.
    % Keeps the same nominal pass/stop edges but is not exactly the paper's
    % equiripple solution.
    b_bp = af_design_bandpass_fir_windowed_sinc(N, fs, f_pass_1, f_pass_2);
    b_bs = af_spectral_invert_fir(b_bp);
    if isempty(af_warned_windowed_fallback) || ~af_warned_windowed_fallback
        warning('AF:NoEquirippleDesigner', ...
                ['firpm/remez not found. Using windowed-sinc FIR fallback (approximate, not exact paper equiripple ', ...
                 'design).']);
        af_warned_windowed_fallback = true;
    end
end


% af_design_visible_variant_filters - Design wide-band visible-artifact FIR pair using A_VIS_* parameters.
%
%   Same six-point frequency-grid structure as af_design_paper_equiripple_filters but
%   uses A_VIS_ENF_HALF_BW_HZ, A_VIS_TRANSITION_BW_HZ, A_VIS_FIR_ORDER, and
%   A_VIS_FIR_WEIGHTS so that attack artifacts are deliberately exaggerated for
%   demonstration and plotting purposes. Falls back to fir1 if firpm fails.
%
%   Inputs:
%     fs        - Sampling rate in Hz.
%     center_hz - Centre frequency of the ENF band (Hz).
%     cfg       - AF_CFG struct with A_VIS_* filter parameters.
%
%   Outputs:
%     b_bp     - Bandpass FIR coefficient vector.
%     b_bs     - Bandstop FIR coefficient vector.
%     is_valid - Logical; false if the requested band falls outside the realizable range.
%
function [b_bp, b_bs, is_valid] = af_design_visible_variant_filters(fs, center_hz, cfg)
    % Same layout as paper filters but with A_VIS_* parameters for
    % exaggerated/visible artifacts in plots.
    nyq = fs / 2;
    half_bw = af_get_cfg_numeric(cfg, 'A_VIS_ENF_HALF_BW_HZ', 5.0);
    trans_bw = af_get_cfg_numeric(cfg, 'A_VIS_TRANSITION_BW_HZ', 1.0);
    N = round(af_get_cfg_numeric(cfg, 'A_VIS_FIR_ORDER', 350));
    W_raw = [1 1 1];
    if isfield(cfg, 'A_VIS_FIR_WEIGHTS') && numel(cfg.A_VIS_FIR_WEIGHTS) == 3
        W_raw = double(cfg.A_VIS_FIR_WEIGHTS(:).');
    end

    f_pass_1 = center_hz - half_bw;
    f_pass_2 = center_hz + half_bw;
    f_stop_1 = f_pass_1 - trans_bw;
    f_stop_2 = f_pass_2 + trans_bw;

    is_valid = (f_stop_1 > 0) && (f_stop_2 < nyq) && ...
               (f_stop_1 < f_pass_1) && (f_pass_2 < f_stop_2);
    if ~is_valid
        b_bp = [];
        b_bs = [];
        return;
    end

    f_edges = [0, f_stop_1, f_pass_1, f_pass_2, f_stop_2, nyq] / nyq;
    a_bp = [0 0 1 1 0 0];
    a_bs = [1 1 0 0 1 1];

    try
        b_bp = firpm(N, f_edges, a_bp, W_raw);
        b_bs = firpm(N, f_edges, a_bs, W_raw);
    catch
        b_bp = fir1(N, [f_pass_1 f_pass_2]/nyq, 'bandpass');
        b_bs = fir1(N, [f_pass_1 f_pass_2]/nyq, 'stop');
    end
end


% af_apply_fir_zero_phase - Apply an FIR filter with zero-phase response (H(f)^2 magnitude).
%
%   Attempts filtfilt for efficient zero-phase filtering; falls back to manual
%   forward-backward application via af_forward_backward_filter if filtfilt is
%   unavailable or the signal is too short (fewer than 3 * filter order samples).
%
%   Inputs:
%     x - Input signal vector.
%     b - FIR filter coefficient vector.
%
%   Outputs:
%     y - Zero-phase filtered signal (same length as x).
%
function y = af_apply_fir_zero_phase(x, b)
    x = x(:);
    if isempty(x) || isempty(b)
        y = x;
        return;
    end

    % Use zero-phase application for aligned before/after plotting.
    try
        if exist('filtfilt', 'file') == 2 && numel(x) > 3 * (numel(b) - 1)
            y = filtfilt(b, 1, x);
        else
            y = af_forward_backward_filter(b, x);
        end
    catch
        y = af_forward_backward_filter(b, x);
    end
    y = y(:);
end


% af_forward_backward_filter - Apply a filter in forward then reverse direction.
%
%   Manual fallback zero-phase implementation used when filtfilt is not available.
%   Applies filter(b,1,x) forward, then filter(b,1,flip(y_forward)) backward.
%   Produces H(f)^2 magnitude response (same as filtfilt for FIR filters).
%
%   Inputs:
%     b - Filter coefficient vector.
%     x - Input signal vector.
%
%   Outputs:
%     y - Zero-phase filtered output (same length as x).
%
function y = af_forward_backward_filter(b, x)
    x = x(:);
    y_f = filter(b, 1, x);
    y_b = filter(b, 1, flipud(y_f));
    y = flipud(y_b);
end


% af_design_bandpass_fir_windowed_sinc - Design a bandpass FIR using windowed-sinc method.
%
%   Last-resort fallback used when neither firpm nor remez is available. Constructs
%   the ideal bandpass impulse response as the difference of two lowpass sincs, then
%   applies a Hamming window and normalises gain at the band centre. Uses only base
%   MATLAB functions (no Signal Processing Toolbox required).
%
%   Inputs:
%     N  - Filter order (even preferred; odd order is incremented by 1).
%     fs - Sampling rate in Hz.
%     f1 - Lower passband edge (Hz).
%     f2 - Upper passband edge (Hz).
%
%   Outputs:
%     b_bp - Bandpass FIR row coefficient vector (1 x N+1).
%
function b_bp = af_design_bandpass_fir_windowed_sinc(N, fs, f1, f2)
    % Base-MATLAB FIR fallback (linear-phase, windowed-sinc).
    % N is filter order, so number of taps is N+1.
    if mod(N, 2) ~= 0
        % Prefer even order for Type-I linear-phase FIR.
        N = N + 1;
    end

    M = N;
    n = (0:M).';
    alpha = M / 2;
    k = n - alpha;

    % Normalize frequencies to cycles/sample (0..0.5)
    fc1 = max(0, f1 / fs);
    fc2 = min(0.5, f2 / fs);
    if ~(fc1 < fc2)
        b_bp = zeros(M+1, 1);
        b_bp(round(alpha)+1) = 1;
        return;
    end

    % Ideal lowpass impulse response helper (cutoff in cycles/sample).
    h_lp_2 = 2 * fc2 * af_sincn(2 * fc2 * k);
    h_lp_1 = 2 * fc1 * af_sincn(2 * fc1 * k);
    h_bp_ideal = h_lp_2 - h_lp_1;

    % Hamming window (implemented locally to avoid toolbox dependency).
    w = 0.54 - 0.46 * cos(2*pi*n / M);
    b_bp = h_bp_ideal .* w;

    % Normalize bandpass by its gain at the band center.
    f0 = (f1 + f2) / 2;
    ejw = exp(-1j * 2*pi * (f0/fs) * (0:M));
    H0 = sum((b_bp(:).') .* ejw);
    g0 = abs(H0);
    if isfinite(g0) && g0 > 0
        b_bp = b_bp / g0;
    end
    b_bp = b_bp(:).';
end


% af_spectral_invert_fir - Convert a bandpass FIR into a complementary bandstop FIR.
%
%   Performs spectral inversion by negating the bandpass coefficients and adding 1
%   to the centre tap, yielding a bandstop filter with the same transition bandwidth
%   and stopband characteristics as the original bandpass.
%
%   Inputs:
%     b_bp - Bandpass FIR coefficient vector.
%
%   Outputs:
%     b_bs - Complementary bandstop FIR coefficient vector (same length).
%
function b_bs = af_spectral_invert_fir(b_bp)
    b_bp = b_bp(:).';
    b_bs = -b_bp;
    mid = floor(numel(b_bs)/2) + 1;
    b_bs(mid) = b_bs(mid) + 1;
end


% af_sincn - Evaluate the normalized sinc function sin(pi*x) / (pi*x).
%
%   Returns 1 where x is zero (limit value). Used by af_design_bandpass_fir_windowed_sinc
%   to construct the ideal lowpass impulse response without requiring sinc() from
%   the Signal Processing Toolbox.
%
%   Inputs:
%     x - Numeric array of any shape.
%
%   Outputs:
%     y - Normalized sinc values; same shape as x.
%
function y = af_sincn(x)
    % Normalized sinc: sin(pi*x)/(pi*x)
    y = ones(size(x));
    idx = abs(x) > eps;
    y(idx) = sin(pi * x(idx)) ./ (pi * x(idx));
end


% af_plot_attack_results - Generate per-attack FFT overlay and spectrogram pair figures.
%
%   For each enabled attack variant (A0_2, A1_2, A2_2), saves three figures:
%   a full-spectrum (0-500 Hz) FFT overlay, a near-ENF (40-80 Hz) FFT overlay,
%   and a side-by-side spectrogram pair. All figures are exported via
%   af_save_figure_svg_and_close. Nested local functions handle FFT computation,
%   single figure assembly, and spectrogram construction to keep the outer body
%   readable. No figure is returned; all output is written to disk.
%
%   Inputs:
%     x_authentic  - Authentic signal vector.
%     x_case_A0_2  - A0 ENF-removed signal (may be empty if runA0_2 is false).
%     x_case_A1_2  - A1 noise-fill signal (may be empty if runA1_2 is false).
%     x_case_A2_2  - A2 synthetic-embedding signal (may be empty if runA2_2 is false).
%     fs           - Sampling rate in Hz.
%     cfg          - AF_CFG struct (needs NOMINAL_FREQ_HZ, A0_ENF_HALF_BW_HZ,
%                   A0_TRANSITION_BW_HZ, and all figure-export fields).
%     runA0_2      - Logical; true to generate A0_2 figures.
%     runA1_2      - Logical; true to generate A1_2 figures.
%     runA2_2      - Logical; true to generate A2_2 figures.
%
%   Outputs:
%     (none) - Figures are saved to disk and closed.
%
function af_plot_attack_results(x_authentic, x_case_A0_2, x_case_A1_2, x_case_A2_2, ...
                                fs, cfg, runA0_2, runA1_2, runA2_2)
    if nargin < 7,  runA0_2 = false; end
    if nargin < 8,  runA1_2 = false; end
    if nargin < 9,  runA2_2 = false; end
    if nargin < 2, x_case_A0_2 = []; end
    if nargin < 3, x_case_A1_2 = []; end
    if nargin < 4, x_case_A2_2 = []; end
    % Generates individual attack plots (no multi-attack subplot dashboard).
    % Harmonic magnitudes for 0-500 Hz views are written to diary log.
    FONT_SZ = 18;
    TITLE_SZ = 20;
    FIG_W = 1400;
    FIG_H = 520;
    SPEC_H = 620;

    function [f_vec, Xdb] = fft_db(x)
        N = length(x);
        X = abs(fft(x(:), N));
        X = X(1:floor(N/2)+1);
        f_vec = (0:floor(N/2)) * fs / N;
        Xdb = 20*log10(X + eps);
    end

    function save_and_close(fig, name)
        if isempty(fig) || ~ishandle(fig)
            return;
        end
        af_save_figure_svg_and_close(fig, cfg, name);
    end

    function fig = plot_fft_overlay_2trace(f, Xdb1, Xdb2, xlims, title_str, label1, label2)
        fig = figure('Position', [100 100 FIG_W FIG_H], 'Color', 'w');
        plot(f, Xdb1, 'b', 'LineWidth', 1.2);
        hold on;
        plot(f, Xdb2, 'r--', 'LineWidth', 1.2);
        hold off;
        xlim(xlims);
        grid on;
        xlabel('Frequency (Hz)', 'FontSize', FONT_SZ);
        ylabel('Magnitude (dB)', 'FontSize', FONT_SZ);
        title(title_str, 'FontSize', TITLE_SZ, 'Interpreter', 'none');
        legend({label1, label2}, 'Location', 'best');
        set(gca, 'FontSize', FONT_SZ);
    end

    function [SdB, F, T] = stft_db(x)
        win_sec = 4.0;
        ov_frac = 0.75;
        nfft_sp = 32768;
        win_samp = round(win_sec * fs);
        hop_samp = round(win_samp * (1 - ov_frac));
        win_samp = max(32, win_samp);
        hop_samp = max(1, hop_samp);
        [S, F, T] = spectrogram(x(:), win_samp, win_samp - hop_samp, nfft_sp, fs);
        SdB = 20*log10(abs(S) + eps);
    end

    function fig = plot_spec_pair(x1, x2, xlims, title1, title2, cax)
        [SdB1, F, T1] = stft_db(x1);
        [SdB2, ~, T2] = stft_db(x2);
        nT = min(numel(T1), numel(T2));
        T = T1(1:nT);
        SdB1 = SdB1(:, 1:nT);
        SdB2 = SdB2(:, 1:nT);
        mask = (F >= xlims(1)) & (F <= xlims(2));
        if ~any(mask)
            fig = [];
            return;
        end

        fig = figure('Position', [100 100 2*FIG_W SPEC_H], 'Color', 'w');
        use_tiled = (exist('tiledlayout', 'file') == 2) || (exist('tiledlayout', 'builtin') == 5);
        if use_tiled
            tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); %#ok<NASGU>
            ax1 = nexttile(1);
        else
            ax1 = subplot(1, 2, 1);
        end
        imagesc(ax1, T, F(mask), SdB1(mask, :));
        axis(ax1, 'xy');
        colormap(ax1, jet);
        xlabel(ax1, 'Time (s)', 'FontSize', FONT_SZ);
        ylabel(ax1, 'Frequency (Hz)', 'FontSize', FONT_SZ);
        title(ax1, title1, 'FontSize', TITLE_SZ, 'Interpreter', 'none');
        set(ax1, 'FontSize', FONT_SZ);
        ylim(ax1, xlims);
        if ~isempty(cax) && numel(cax) == 2 && all(isfinite(cax)) && cax(2) > cax(1)
            caxis(ax1, cax);
        end

        if use_tiled
            ax2 = nexttile(2);
        else
            ax2 = subplot(1, 2, 2);
        end
        imagesc(ax2, T, F(mask), SdB2(mask, :));
        axis(ax2, 'xy');
        colormap(ax2, jet);
        xlabel(ax2, 'Time (s)', 'FontSize', FONT_SZ);
        ylabel(ax2, 'Frequency (Hz)', 'FontSize', FONT_SZ);
        title(ax2, title2, 'FontSize', TITLE_SZ, 'Interpreter', 'none');
        set(ax2, 'FontSize', FONT_SZ);
        ylim(ax2, xlims);
        if ~isempty(cax) && numel(cax) == 2 && all(isfinite(cax)) && cax(2) > cax(1)
            caxis(ax2, cax);
        end
        colorbar(ax2);
    end

    [f_h, Xdb_h] = fft_db(x_authentic);

    % Shared caxis from the authentic trace for variant spectrogram comparisons.
    cax_a0_shared_40_80 = [];
    try
        [S_h, F_h, ~] = spectrogram(x_authentic(:), round(4.0*fs), round(4.0*fs*0.75), 32768, fs);
        SdB_h = 20*log10(abs(S_h) + eps);
        % Compute 40-80 Hz caxis from off-ENF background only.
        enf_excl_lo = cfg.NOMINAL_FREQ_HZ - (cfg.A0_ENF_HALF_BW_HZ + cfg.A0_TRANSITION_BW_HZ) - 2;
        enf_excl_hi = cfg.NOMINAL_FREQ_HZ + (cfg.A0_ENF_HALF_BW_HZ + cfg.A0_TRANSITION_BW_HZ) + 2;
        mask_bg = (F_h >= 40) & (F_h <= 80) & ~((F_h >= enf_excl_lo) & (F_h <= enf_excl_hi));
        if any(mask_bg)
            v_bg = SdB_h(mask_bg, :);
            bg_lo = prctile(v_bg(:), 2);
            bg_hi = prctile(v_bg(:), 98);
            cax_a0_shared_40_80 = [bg_lo - 5, bg_hi + 20];
        end
    catch
        % Keep defaults (empty) and let plots auto-scale if the authentic STFT fails.
    end

    if runA0_2 && ~isempty(x_case_A0_2)
        [~, Xdb_a0_2] = fft_db(x_case_A0_2);
        save_and_close(plot_fft_overlay_2trace(f_h, Xdb_h, Xdb_a0_2, [0 500], 'A0_2 Overlay FFT 0-500Hz', 'Authentic', 'A0_2'), 'A0_2_overlay_fft_0_500hz');
        save_and_close(plot_fft_overlay_2trace(f_h, Xdb_h, Xdb_a0_2, [40 80], 'A0_2 Overlay FFT 40-80Hz', 'Authentic', 'A0_2'), 'A0_2_overlay_fft_40_80hz');
        save_and_close(plot_spec_pair(x_authentic, x_case_A0_2, [40 80], 'Authentic 40-80Hz', 'A0_2 40-80Hz', cax_a0_shared_40_80), 'A0_2_pair_spec_40_80hz');
    end

    if runA1_2 && ~isempty(x_case_A1_2)
        [~, Xdb_a1_2] = fft_db(x_case_A1_2);
        save_and_close(plot_fft_overlay_2trace(f_h, Xdb_h, Xdb_a1_2, [0 500], 'A1_2 Overlay FFT 0-500Hz', 'Authentic', 'A1_2'), 'A1_2_overlay_fft_0_500hz');
        save_and_close(plot_fft_overlay_2trace(f_h, Xdb_h, Xdb_a1_2, [40 80], 'A1_2 Overlay FFT 40-80Hz', 'Authentic', 'A1_2'), 'A1_2_overlay_fft_40_80hz');
        save_and_close(plot_spec_pair(x_authentic, x_case_A1_2, [40 80], 'Authentic 40-80Hz', 'A1_2 40-80Hz', cax_a0_shared_40_80), 'A1_2_pair_spec_40_80hz');
    end

    if runA2_2 && ~isempty(x_case_A2_2)
        [~, Xdb_a2_2] = fft_db(x_case_A2_2);
        save_and_close(plot_fft_overlay_2trace(f_h, Xdb_h, Xdb_a2_2, [0 500], 'A2_2 Overlay FFT 0-500Hz', 'Authentic', 'A2_2'), 'A2_2_overlay_fft_0_500hz');
        save_and_close(plot_fft_overlay_2trace(f_h, Xdb_h, Xdb_a2_2, [40 80], 'A2_2 Overlay FFT 40-80Hz', 'Authentic', 'A2_2'), 'A2_2_overlay_fft_40_80hz');
        save_and_close(plot_spec_pair(x_authentic, x_case_A2_2, [40 80], 'Authentic 40-80Hz', 'A2_2 40-80Hz', cax_a0_shared_40_80), 'A2_2_pair_spec_40_80hz');
    end

end


function s = af_cfg_get_string_or_default(cfg, field_name, default_val)
    s = char(default_val);
    if nargin < 3 || isempty(default_val)
        s = '';
    end
    if nargin < 1 || isempty(cfg) || ~isstruct(cfg)
        return;
    end
    if ~isfield(cfg, field_name)
        return;
    end
    v = cfg.(field_name);
    try
        if isstring(v) && isscalar(v) && strlength(v) > 0
            s = char(v);
        elseif ischar(v) && ~isempty(v)
            s = v;
        end
    catch
        % keep default
    end
end


% af_fft_one_sided_no_dc_core - Compute the one-sided complex FFT spectrum with DC bin removed.
%
%   Computes the Nfft-point FFT of x, retains the positive-frequency half
%   (bins 0..Nfft/2), and removes the DC bin (f=0) before returning.
%   Used as the shared low-level FFT primitive by af_fft_mag_abs_no_dc and
%   related helpers.
%
%   Inputs:
%     x    - Signal vector (converted to double internally).
%     fs   - Sampling rate in Hz (used to build the frequency axis).
%     Nfft - FFT size (must be >= 2).
%
%   Outputs:
%     f_fft - Frequency axis vector (Hz), DC bin excluded.
%     X_pos - Complex FFT coefficients corresponding to f_fft.
%
function [f_fft, X_pos] = af_fft_one_sided_no_dc_core(x, fs, Nfft)
    f_fft = [];
    X_pos = [];
    x = double(x(:));
    if numel(x) < 2 || ~isfinite(fs) || fs <= 0 || ~isfinite(Nfft) || Nfft < 2
        return;
    end

    X = fft(x, Nfft);
    X = X(1:floor(Nfft/2)+1);
    f_all = (0:floor(Nfft/2)).' * (fs / Nfft);

    idx = f_all > 0; % remove DC
    f_fft = f_all(idx);
    X_pos = X(idx);
end


function val = af_summarize_values(v, stat_mode)
    v = v(isfinite(v));
    if isempty(v)
        val = NaN;
        return;
    end
    switch lower(stat_mode)
        case 'mean'
            val = mean(v);
        otherwise
            val = median(v);
    end
end


% af_fft_mag_abs_no_dc - Compute one-sided FFT magnitude with DC removed and next-power-of-2 zero padding.
%
%   Removes DC bias via detrend, pads to the next power of 2, and delegates to
%   af_fft_one_sided_no_dc_core to obtain the positive-frequency complex spectrum.
%   Returns the absolute (linear) magnitude. Used throughout the script as the
%   standard FFT utility for spectral measurement and reporting.
%
%   Inputs:
%     x  - Signal vector.
%     fs - Sampling rate in Hz.
%
%   Outputs:
%     f_fft   - One-sided frequency axis (Hz), DC excluded.
%     mag_abs - Linear magnitude vector corresponding to f_fft.
%
function [f_fft, mag_abs] = af_fft_mag_abs_no_dc(x, fs)
    x = x(:);
    N = numel(x);
    if N < 2
        f_fft = [];
        mag_abs = [];
        return;
    end

    x = detrend(x, 0);
    Nfft = 2^nextpow2(N);
    [f_fft, X_pos] = af_fft_one_sided_no_dc_core(x, fs, Nfft);
    mag_abs = abs(X_pos);
end


% af_get_aligned_enf_for_plot - Return time-aligned ENF estimate pairs for overlay plotting.
%
%   Builds a per-frame time axis for each ENF vector using hop_sec = (frame_size -
%   overlap_size) / fs, shifts enf2 by lag_sec, finds the overlapping time interval,
%   and interpolates enf2 onto enf1's time grid. Removes frames where either
%   estimate is non-finite. Returns empty arrays if the overlap is shorter than
%   min_overlap_sec.
%
%   Inputs:
%     enf1            - ENF estimate vector for trace 1 (frames x 1, Hz).
%     enf2            - ENF estimate vector for trace 2 (frames x 1, Hz).
%     fs              - Sampling rate of the original signal in Hz.
%     frame_size      - STFT window length in samples.
%     overlap_size    - STFT overlap in samples.
%     lag_sec         - Time shift applied to enf2 time axis (seconds).
%     min_overlap_sec - (optional) Minimum required overlap in seconds (default 0).
%
%   Outputs:
%     t  - Common time vector (seconds) for the overlapping region.
%     y1 - enf1 values at times t.
%     y2 - enf2 values interpolated to times t.
%
function [t, y1, y2] = af_get_aligned_enf_for_plot(enf1, enf2, fs, frame_size, overlap_size, lag_sec, min_overlap_sec)
% enf1/enf2  : ENF estimate vectors (frames x 1)
% lag_sec    : time offset applied to enf2 (positive = enf2 shifted later)
% min_overlap_sec : minimum required overlap; returns empty if not met.
    if nargin < 7
        min_overlap_sec = 0;
    end
    hop_sec = (frame_size - overlap_size) / fs;
    enf1 = enf1(:);
    enf2 = enf2(:);
    t1 = (0:numel(enf1)-1).' * hop_sec;
    t2 = (0:numel(enf2)-1).' * hop_sec + lag_sec;
    t_start = max(t1(1),  t2(1));
    t_end   = min(t1(end), t2(end));
    if (t_end - t_start) < min_overlap_sec
        t = []; y1 = []; y2 = [];
        return;
    end
    idx1 = (t1 >= t_start) & (t1 <= t_end);
    t  = t1(idx1);
    y1 = enf1(idx1);
    y2 = interp1(t2, enf2, t, 'linear', NaN);
    good = isfinite(y1) & isfinite(y2);
    t  = t(good);
    y1 = y1(good);
    y2 = y2(good);
end
