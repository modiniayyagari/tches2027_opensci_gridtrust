%% ENF Signature Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% This script implements the signal processing pipeline for the application
% "Geolocation for Semiconductor Chips using Electric Network Frequency Signatures".
%
% The primary goal is to validate the feasibility of extracting Electric
% Network Frequency (ENF) signatures from DC-powered hardware by comparing
% a sensed trace against a ground-truth reference.
%
% Analysis Steps:
%   1.  INPUTS: The script takes two time-synchronized inputs (both must be .wav):
%       - A 'ground-truth reference trace' captured from the AC mains. (.wav file format)
%       - A 'sensed trace' captured from the experimental FPGA board (.wav file format).
%         This can be either an ambient EM trace or a board power trace.
%
%   2.  SPECTROGRAM GENERATION: The Short-Time Fourier Transform (STFT) is
%       applied to both traces to generate high-resolution spectrograms,
%       visualizing the harmonic content over time.
%
%   3.  ENF ESTIMATION: A weighted average PMF is used to extract the instantaneous ENF
%       signature from each spectrogram.
%
%   4.  OUTPUTS & CORRELATION: Finally, the script compares the ENF signature
%       from the sensed trace against the signature from the ground-truth
%       reference. It calculates the Pearson correlation coefficient and
%       generates the final temporally aligned plots to visually and quantitatively
%       assess the match.
%
% Dependencies:
%   proc_enf_analysis - Pre-compiled ENF extraction and correlation function.
%                       Must be on the MATLAB path before running this script.
%
% Usage:
%   Run from the exp_scripts/ directory. The artifact data root is resolved
%   automatically relative to this script file's location.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENF Analysis Configuration Flags:
% The parameters in this section control what is analyzed, how the STFT spectrogram
% is computed, how the ENF frequency is estimated, and how results are displayed.
% Update the four variables in the 'Input trace file information' block for each
% new run; all other parameters can typically remain at their default values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chunk-processing configuration
%
% SCRP_CFG_DO_CHUNK_PROCESSING:
%   true:  Split the full .wav traces into smaller segments and analyze them.
%   false: Skip segmenting and analyze only the full-length .wav traces.
SCRP_CFG_DO_CHUNK_PROCESSING = false;
%
% SCRP_CFG_ENABLE_CHUNK_PLOTS:
%   true:  Generate plots for each segment pair.
%   false: Compute segment PCCs without making per-segment plots.
SCRP_CFG_ENABLE_CHUNK_PLOTS = false;
%
% SCRP_CFG_CHUNK_PLOT_MIN_OVERLAP_SEC:
%   Used only if SCRP_CFG_ENABLE_CHUNK_PLOTS=true.
SCRP_CFG_CHUNK_PLOT_MIN_OVERLAP_SEC = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input trace file information
% This script is designed to analyze ONE pair of traces at a time. Manually
% update the variables below for EACH run. If both traces are in the same
% leaf folder, set file_1_path and file_2_path to the same value.
%
% 1. `file_1_path`: Relative path ending in "/" for the AC mains reference.
% 2. `file_2_path`: Relative path ending in "/" for the sensed trace.
% 3. `file_1_name`: AC mains reference file name without extension.
% 4. `file_2_name`: Sensed trace file name without extension.
%
% --- EXAMPLES ---
%
% % Example 1: Analyzing a 'TREND' experiment trace (WK01/WEND/SUN/EVEN/T01)
% file_1_path = "../exp_inputs/TREND/WK01/WEND/SUN/EVEN/T01/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 2: Analyzing an FPGA ambient EM trace (FPGA/EM_TRACES/CW305/DC_1)
% file_1_path = "../exp_inputs/FPGA/EM_TRACES/CW305/DC_1/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 2a: Analyzing an FPGA ambient EM trace for ON-OFF-ON EXP(FPGA/EM_TRACES/SAKU/ON_OFF_ON)
% file_1_path = "../exp_inputs/FPGA/EM_TRACES/SAKU/ON_OFF_ON/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 3: Analyzing an FPGA board power trace (FPGA/POW_TRACES/CW305/DC_1)
% file_1_path = "../exp_inputs/FPGA/POW_TRACES/CW305/DC_1/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_pow_trace_dc_a";  % or "fpga_pow_trace_dc_b"
%
% % Example 4: Analyzing a 'MULTI' (US_60) experiment trace (MULTI/US_60/AUG/WED/T02)
% % File names differ from the standard naming convention in this folder.
% file_1_path = "../exp_inputs/MULTI/US_60/AUG/WED/T02/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac_egrid_citya_lab";
% file_2_name = "fpga_em_trace_dc_egrid_citya_lab";
%
% % Example 5: Analyzing a 'MULTI' (DE_50) experiment trace (MULTI/DE_50/AUG/WED/T01)
% file_1_path = "../exp_inputs/MULTI/DE_50/AUG/WED/T01/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac_citya_lab";
% file_2_name = "fpga_em_trace_dc_citya_lab";
%
% % Example 6: Analyzing an 'SRV_L' experiment trace (SRV_L/SPSU/WED/T05)
% file_1_path = "../exp_inputs/SRV_L/SPSU/WED/T05/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";

% Set the variables below to point to the trace pair you want to analyze.
% file_1_path = "../exp_inputs/FPGA/EM_TRACES/CW305/DC_1/";
% file_2_path = file_1_path;
file_1_path = "../exp_inputs/FPGA/EM_TRACES/SAKU/DC_1/";
file_2_path = file_1_path;
file_1_name = "mains_pow_trace_ac";
file_2_name = "fpga_em_trace_dc";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram Computation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental power grid frequency settings
nominal_freq_arr      = [50 60];

% --- Multi-Location Experiment Note (DE_50 / US_60) ---
% When analyzing traces from the 'MULTI' directory, you MUST manually
% update the `nominal_freq_*` variables below to match the fundamental
% frequency of the grid from where the trace was recorded.
% - For US_60 (United States) traces: set nominal_freq = nominal_freq_arr(2); (60 Hz)
% - For DE_50 (Germany) traces:     set nominal_freq = nominal_freq_arr(1); (50 Hz)

%For reference trace
nominal_freq_1        = nominal_freq_arr(2);
harmonics_arr_1       = (1:7)*nominal_freq_1;

%For sensed trace
nominal_freq_2        = nominal_freq_arr(2);
harmonics_arr_2       = (1:7)*nominal_freq_2;

% STFT compute param settings
frame_size_arr      = (1:12)*1000;
frame_size          = frame_size_arr(8);                %8000ms window
nfft_arr            = 2.^(10:20);
nfft                = nfft_arr(6);                      %2^15 = 32768 pts
overlap_size_arr    = 0:0.1:0.9;
overlap_size        = overlap_size_arr(1)*frame_size;   %non-overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Estimation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency to be estimated
trace_1_est_freq = harmonics_arr_1(1);                    %1st harmonic= 60Hz
trace_2_est_freq = harmonics_arr_2(1);                    %1st harmonic= 60Hz

%Frequency estimation method options:
%1: weighted average (pmf power = 3)
%2: spectrum combining (quad interp)
trace_1_freq_est_method = 1;   %For reference trace
trace_1_freq_est_spec_comb_harmonics = [60 120 180 240 300 360 420];

trace_2_freq_est_method = 1;   %For sensed trace
trace_2_freq_est_spec_comb_harmonics = [60 120 180 240 300 360 420];

tempo_align_corr_vals = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Plotting Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The colormap for the figures is set to jet.
set(0,'DefaultFigureColormap', jet)

%Titles for MATLAB plots
trace_1_plot_title = "Reference AC Mains Power Trace";
trace_2_plot_title = "Sensed FPGA Trace";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% File Path Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the full .wav file paths for the reference and sensed traces by
% concatenating each configured leaf-folder path, file base name, and the
% '.wav' extension. Both inputs must already be in .wav format; no format
% conversion is performed here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_file_1_path_wav = file_1_path + file_1_name + ".wav";
full_file_2_path_wav = file_2_path + file_2_name + ".wav";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Chunk Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optionally split the long traces into custom time segments and analyze
% each reference/sensed pair separately.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SCRP_CFG_DO_CHUNK_PROCESSING
    % Custom segment boundaries in seconds:
    % [0 50 120 Inf] creates:
    %   segment 1: 0   to 50 seconds
    %   segment 2: 50  to 120 seconds
    %   segment 3: 120 to end of file
    segment_edges_secs = [0 50 120 Inf];

    output_filename_w_chunk_1 = file_1_path + file_1_name + "_tr";
    fprintf('\nINFO: For trace 1 splitting "%s" into custom time segments...\n', full_file_1_path_wav);
    trace_1_created_files = proc_split_trace_segments( ...
        full_file_1_path_wav, segment_edges_secs, output_filename_w_chunk_1);

    output_filename_w_chunk_2 = file_2_path + file_2_name + "_tr";
    fprintf('\nINFO: For trace 2 splitting "%s" into custom time segments...\n', full_file_2_path_wav);
    trace_2_created_files = proc_split_trace_segments( ...
        full_file_2_path_wav, segment_edges_secs, output_filename_w_chunk_2);
else
    disp('INFO: Skipping trace chunking as per SCRP_CFG_DO_CHUNK_PROCESSING flag.');
    trace_1_created_files = {};
    trace_2_created_files = {};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Main Analysis: ENF Extraction and Correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the full ENF analysis pipeline on the configured reference and sensed
% trace pair. All signal processing is handled by proc_enf_analysis, which:
%   - Computes the STFT spectrogram for both traces
%   - Extracts the instantaneous ENF time series from each spectrogram
%   - Temporally aligns the two ENF traces and computes Pearson correlation
%   - Generates all intermediate spectrograms and the final correlation plot
% The returned value is the best (temporally aligned) Pearson correlation (r).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nINFO: ENF Analysis on complete reference and sensed traces ...\n');
matched_corr_val_full = proc_enf_analysis(full_file_1_path_wav, full_file_2_path_wav, ...
                                          nfft, frame_size, overlap_size, ...
                                          harmonics_arr_1, nominal_freq_1, harmonics_arr_2, nominal_freq_2,...
                                          trace_1_freq_est_method, trace_1_est_freq, trace_1_freq_est_spec_comb_harmonics, ...
                                          trace_2_freq_est_method, trace_2_est_freq, trace_2_freq_est_spec_comb_harmonics, ...
                                          trace_1_plot_title, trace_2_plot_title, true);

tempo_align_corr_vals{end+1} = matched_corr_val_full;

if SCRP_CFG_DO_CHUNK_PROCESSING
    if numel(trace_1_created_files) == numel(trace_2_created_files)
        fprintf('\nINFO: ENF Analysis (pairwise) on %d trace chunks for reference and sensed traces ...\n', numel(trace_1_created_files));
        for i = 1:numel(trace_1_created_files)
            file1 = trace_1_created_files{i};
            file2 = trace_2_created_files{i};
            fprintf('\nINFO: Processing pair %d: "%s" and "%s"\n', i, file1, file2);
            matched_corr_val_chunk = proc_enf_analysis(file1, file2, ...
                                                       nfft, frame_size, overlap_size, ...
                                                       harmonics_arr_1, nominal_freq_1, harmonics_arr_2, nominal_freq_2,...
                                                       trace_1_freq_est_method, trace_1_est_freq, trace_1_freq_est_spec_comb_harmonics, ...
                                                       trace_2_freq_est_method, trace_2_est_freq, trace_2_freq_est_spec_comb_harmonics, ...
                                                       trace_1_plot_title, trace_2_plot_title, SCRP_CFG_ENABLE_CHUNK_PLOTS, ...
                                                       'MinOverlapSec', SCRP_CFG_CHUNK_PLOT_MIN_OVERLAP_SEC);
            tempo_align_corr_vals{end+1} = matched_corr_val_chunk;
        end

        fprintf('\nINFO: Pearsons Correlation Coefficient for the temporally aligned reference and sensed traces (complete trace): %.8f \n', cell2mat(tempo_align_corr_vals(1)));
        fprintf('\nINFO: Pearsons Correlation Coefficients for the temporally aligned reference and sensed traces (chunk-wise):\n');
        disp(tempo_align_corr_vals(2:end));
        mean_corr = mean(cell2mat(tempo_align_corr_vals(2:end)));
        fprintf('INFO: Mean Pearsons Correlation Coefficient for the temporally aligned reference and sensed traces (all chunks): %.8f\n', mean_corr);

        fprintf('\nINFO: ENF Analysis Post-Processing complete.\n');
    else
        error('The number of files in trace_1_created_files and trace_2_created_files must be the same.');
    end
else
    fprintf('\nINFO: Pearsons Correlation Coefficient for the temporally aligned reference and sensed traces (complete trace): %.8f \n', matched_corr_val_full);
    fprintf('\nINFO: ENF Analysis complete.\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Local Helper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split a WAV file using user-defined time boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function created_files = proc_split_trace_segments(inputFile, segment_edges_secs, output_filename_suffix)
    if nargin < 3
        output_filename_suffix = 'tr';
    end

    created_files = {};

    try
        [y, Fs] = audioread(inputFile);
        [totalSamples, ~] = size(y);
        totalDurationSecs = totalSamples / Fs;

        segment_edges_secs(~isfinite(segment_edges_secs)) = totalDurationSecs;

        if segment_edges_secs(1) < 0
            error('Segment start time cannot be negative.');
        end

        if any(diff(segment_edges_secs) <= 0)
            error('segment_edges_secs must be strictly increasing.');
        end

        if segment_edges_secs(1) ~= 0
            warning('First segment does not start at 0 seconds.');
        end

        segment_edges_secs(end) = min(segment_edges_secs(end), totalDurationSecs);
        numSegments = numel(segment_edges_secs) - 1;

        for i = 1:numSegments
            startSec = segment_edges_secs(i);
            endSec   = segment_edges_secs(i + 1);

            startSample = floor(startSec * Fs) + 1;
            endSample   = min(floor(endSec * Fs), totalSamples);

            if endSample <= startSample
                warning('Skipping segment %d because it has no samples.', i);
                continue;
            end

            currentSegmentData = y(startSample:endSample, :);
            outputFileName = sprintf('%s_%d.wav', output_filename_suffix, i);
            audiowrite(outputFileName, currentSegmentData, Fs);

            fprintf('INFO: Created segment %d: %.2f to %.2f sec -> %s\n', ...
                    i, startSec, endSec, outputFileName);

            created_files{end+1} = outputFileName;
        end

    catch ME
        fprintf('ERROR: Could not process the file "%s".\n', inputFile);
        fprintf('MATLAB error message: %s\n', ME.message);
        rethrow(ME);
    end
end
