# TCHES 2027 Cycle 1 Submission: “Dude, Where's My Chip? Electric Network Frequency Signatures for Chip Geolocation”

## 1. Introduction

This repository contains the artifacts for the paper titled **“Dude, Where's My Chip? Electric Network Frequency Signatures for Chip Geolocation,”** submitted to **TCHES 2027 Cycle 1**. The artifacts provided here allow for independent verification of our main experimental results on extracting and using Electric Network Frequency (ENF) signatures from the local environment of **DC-powered PCBs / FPGA boards**.

The repository includes:

- The raw ENF datasets used in the paper:
  - Baseline ambient EM feasibility measurements (Experiment 1)
  - FPGA-based ambient EM traces (Experiment 2)
  - FPGA ambient EM ON-OFF-ON control trace (Experiment 2)
  - FPGA board power-trace measurements (Experiment 2)
  - Temporal reliability measurements at a single grid location (Experiment 3).
  - Server-room robustness measurements (Experiment 4).
  - Multi-location cross-grid validation measurements (Experiment 5).
- MATLAB scripts that implement the ENF extraction, correlation analysis, statistical summaries, and attack analyses.
- Attack-analysis scripts for stale replay, relay-delay sensitivity, and synthetic ENF forgery attack (Security Analysis Section)
- A compiled MATLAB analysis core (`proc_enf_analysis.p`) used by the top-level scripts.
- The additional experiment companion PDF `enf_based_chip_geolocation_additional_experiment_artifacts.pdf`, which provides the supplementary visual artifacts for the added control experiment material.

*Note:* Our core analysis function, originally developed under authors' university policy that restricts source release, is provided here as a compiled MATLAB function (`proc_enf_analysis.p`) rather than in source form. The top-level scripts in this repository call that compiled function so that all analyses in the paper can be reproduced without exposing proprietary internals.


## 2. Folder and File Hierarchy

At a high level, the repository is organized as follows:

- **`exp_inputs/`** – All ENF input traces (ambient EM + board power traces + mains references) used in the experiments.

  - **`BASE/`** – Baseline ambient EM feasibility experiment (no FPGA workload / simple setup).
    - **`NO_FB/`**
      - `fpga_em_trace_dc.wav` – Ambient EM trace captured with only picoscope.
      - `mains_pow_trace_ac.wav` – Ground-truth mains reference trace.
    - **`W_FB/`**
      - Same file naming as `NO_FB`, but with the setup placed inside a shielding (e.g., Faraday bag) to test isolation effects.

  - **`FPGA/`** – FPGA-board ENF feasibility experiment.
    This folder is split by sensed-trace type so that ambient EM traces and board power traces can be analyzed independently.

    - **`EM_TRACES/`** – FPGA-board ambient EM traces.
      - **`CW305/`**
        - **`DC_1/`–`DC_5/`**
          - `fpga_em_trace_dc.wav` – Ambient EM trace near the CW305 FPGA board.
          - `mains_pow_trace_ac.wav` – Co-recorded mains reference.
      - **`SAKU/`**
        - **`DC_1/`–`DC_5/`**
          - `fpga_em_trace_dc.wav` – Ambient EM near the Sakura-G (or equivalent) FPGA board.
          - `mains_pow_trace_ac.wav` – Co-recorded mains reference.
        - **`ON_OFF_ON/`**
          - `fpga_em_trace_dc.wav` – Ambient EM trace for the Sakura-G ON-OFF-ON control experiment.
          - `mains_pow_trace_ac.wav` – Co-recorded mains reference.
          - Single 200 s control recording with FPGA supply connected from 0-50 s, disconnected from 50-120 s, and reconnected from 120-200 s.

    - **`POW_TRACES/`** – FPGA-board power-trace measurements.
      - **`CW305/`**
        - **`DC_1/`–`DC_5/`**
          - `fpga_pow_trace_dc_a.wav` – FPGA power-rail trace from VIN pin of voltage regulator
          - `fpga_pow_trace_dc_b.wav` – FPGA power-rail trace from VOUT pin of voltage regulator
          - `mains_pow_trace_ac.wav` – Co-recorded mains reference.
      - **`SAKU/`**
        - **`DC_1/`–`DC_5/`**
          - `fpga_pow_trace_dc_a.wav` – FPGA power-rail trace from VIN pin of voltage regulator
          - `fpga_pow_trace_dc_b.wav` – FPGA power-rail trace from VOUT pin of voltage regulator
          - `mains_pow_trace_ac.wav` – Co-recorded mains reference.

  - **`MULTI/`** – **Experiment 5 – ENF Multi-Location Validation.**
    ENF measurements comparing a single local EM trace against mains references from multiple grid regions. The same dataset also supports relay-delay attack analysis.

    - **`US_60/`** – 60 Hz (US) multi-location dataset.
      - **`AUG/`**, **`OCT/`** – Months of collection.
        - `MON/`, `WED/`, `TUE/`, etc. – Calendar day.
          - `T01/`, `T02/`, `T03/` – Measurement folders (each one 200-s trace pair set).
            - `fpga_em_trace_dc_egrid_citya_lab.wav` – Local ambient EM trace near the FPGA.
            - `mains_pow_trace_ac_egrid_citya_lab.wav` – Local mains reference (same grid as EM sensor).
            - `mains_pow_trace_ac_egrid_citya_home.wav` – Additional local mains (home) reference.
            - `mains_pow_trace_ac_egrid_worcester.wav` – US Eastern grid reference.
            - `mains_pow_trace_ac_tgrid_richardson.wav` – US Texas grid reference.
            - `mains_pow_trace_ac_wgrid_tucson.wav` – US Western grid reference.

    - **`DE_50/`** – 50 Hz (European) multi-location dataset.
      - **`AUG/`, `OCT/`** – Months of collection.
        - `WED/`, `THU/`, `TUE/`, … – Calendar day.
          - `T01/`, `T02/`, `T03/` – Measurement folders.
            - `fpga_em_trace_dc_citya_lab.wav` – Local 50 Hz EM trace near FPGA.
            - `mains_pow_trace_ac_citya_lab.wav` – Local 50 Hz mains reference.
            - `mains_pow_trace_ac_dresden.wav` – Remote European mains reference.

    Each `Txx` folder thus contains one EM trace and multiple mains references from different grid regions, enabling the cross-grid correlation matrices used in Experiment 5.

  - **`TREND/`** – **Experiment 3 – ENF Temporal Reliability.**
    Designed to measure how stable ENF correlation is over **Week**, **Day-of-Week**, and **Time-of-Day** at a single grid location. This dataset is also used for stale replay attack analysis.

    - **`WK01/`, `WK02/`** – Two consecutive weeks.
      - **`WDAY/`**, **`WEND/`** – Weekday vs weekend partition.
        - **`WED/`, `THUR/`, `SAT/`, `SUN/`** – Specific days.
          - **`EMRN/`, `MORN/`, `AFTN/`, `EVEN/`** – Time-of-day slots:
            - EMRN ≈ early morning (4 AM).
            - MORN ≈ morning (9 AM).
            - AFTN ≈ afternoon (2 PM).
            - EVEN ≈ evening (7 PM).
          - **`T01/`–`T05/`** – Five repeated trials per (Week, Day, Time-of-Day) condition.
            - `fpga_em_trace_dc.wav` – Ambient EM trace near the FPGA board.
            - `mains_pow_trace_ac.wav` – Mains reference trace.

    This design yields 2 weeks × 4 days × 4 times × 5 trials = 160 measurements.

  - **`SRV_L/`** – **Experiment 4 – ENF Server-Room Robustness.**
    ENF measurements in a high-density computing environment, studying robustness to sensor placement on the server PSU and day-of-week.

    - **`SBOX/`**, **`SPSU/`** – Two sensor locations:
      - For each site:
        - **`MON/`, `WED/`** – Two days.
          - **`T01/`–`T05/`** – Five repeated trials per (Site, Day) condition.
            - `fpga_em_trace_dc.wav` – Ambient EM trace near server PSU/BOX.
            - `mains_pow_trace_ac.wav` – Co-recorded mains reference.

- **`exp_scripts/`** – MATLAB scripts required to run the ENF extraction, correlation, statistical analysis, and attack analysis.

  - `enf_analysis_top_pair_wise.m`
    – Top-level script for basic pair-wise ENF signature analysis (for singular analysis).

  - `enf_analysis_top_reliability_stat.m`
    – Top-level script for **Experiment 3 – Temporal Reliability** (TREND), including stale replay attack analysis.

  - `enf_analysis_top_reliability_stat_server.m`
    – Top-level script for **Experiment 4 – Server-Room Robustness** (SRV_L).

  - `enf_analysis_top_multi_loc_stat.m`
    – Top-level script for **Experiment 5 – Multi-Location Validation** (MULTI), including the relay-delay attack mode.

  - `enf_analysis_top_pair_wise_af.m`
    – Top-level script for pair-wise synthetic forgery attack.

  - `proc_enf_analysis.p`
    – Compiled MATLAB function that implements the core ENF extraction and correlation pipeline shared across all scripts.

- **`README.md`** – This file.


## 3. Experiment Scripts Overview

This section summarizes what each top-level MATLAB script does and how it relates to the paper’s experiments.

### 3.1 Pair-Wise ENF Signature Analysis (Feasibility Pipeline)

**Script:** `enf_analysis_top_pair_wise.m`
**Data roots:** `exp_inputs/*`

**Goal:** Validate the feasibility of extracting ENF from **DC-powered hardware** by comparing a sensed ambient EM trace or board power trace against a mains reference.

**Overview:**

1. **Inputs**
   - `fpga_em_trace_dc*.wav`: Ambient EM trace captured near the FPGA/PCB.
   - `fpga_pow_trace_dc_*.wav`: Board power-trace channel captured from the FPGA setup.
   - `mains_pow_trace_ac*.wav`: Co-recorded mains reference trace.
   - The script uses separate `file_1_path` and `file_2_path` variables:
     - `file_1_path` points to the mains reference folder.
     - `file_2_path` points to the sensed-trace folder.
     - For most datasets these two paths are identical.

   FPGA examples:
   - Ambient EM trace:
     - `file_1_path = "../exp_inputs/FPGA/EM_TRACES/CW305/DC_1/";`
     - `file_2_path = file_1_path;`
     - `file_1_name = "mains_pow_trace_ac";`
     - `file_2_name = "fpga_em_trace_dc";`
   - Ambient EM ON-OFF-ON control trace:
     - `file_1_path = "../exp_inputs/FPGA/EM_TRACES/SAKU/ON_OFF_ON/";`
     - `file_2_path = file_1_path;`
     - `file_1_name = "mains_pow_trace_ac";`
     - `file_2_name = "fpga_em_trace_dc";`
   - Board power trace:
     - `file_1_path = "../exp_inputs/FPGA/POW_TRACES/CW305/DC_1/";`
     - `file_2_path = file_1_path;`
     - `file_1_name = "mains_pow_trace_ac";`
     - `file_2_name = "fpga_pow_trace_dc_a";`
     - Use `fpga_pow_trace_dc_b` to analyze the second power-trace channel.

2. **Spectrogram Generation**
   - Applies Short-Time Fourier Transform (STFT) to both traces to obtain time–frequency spectrograms around the nominal grid frequency and its harmonics.

3. **ENF Estimation**
   - Uses a weighted-harmonic ENF estimator to extract the instantaneous ENF trajectory from each spectrogram.

4. **Correlation & Outputs**
   - Computes the **Pearson correlation coefficient** between the sensed ENF and the mains ENF.
   - Optionally produces aligned ENF plots and summary correlation statistics for the baseline and FPGA experiments.

This script is the simplest entry point if you want to understand the core ENF extraction and correlation pipeline on a single pair of traces.

#### Experiment 2 ON-OFF-ON Control Trace

This additional controlled dataset supports **Experiment 2** by isolating the effect of the energized DC-powered FPGA path on ambient EM ENF sensing. The measurement setup uses the same Sakura-G ambient EM sensing arrangement as the main FPGA ambient EM experiment: an unterminated BNC cable is placed near the Sakura-G FPGA board to capture local ambient EM, while a synchronized AC mains trace is recorded as the reference.

The control recording is stored under `exp_inputs/FPGA/EM_TRACES/SAKU/ON_OFF_ON/` and contains one 200 s trace pair:

- 0-50 s: FPGA supply connected.
- 50-120 s: FPGA supply physically disconnected from the AC mains source.
- 120-200 s: FPGA supply reconnected.

The sensing probe and acquisition setup remain unchanged during the full recording. When analyzed with the pair-wise script, the two ON intervals show strong ENF agreement with the AC mains reference, while the middle OFF interval drops toward the baseline ambient EM behavior from Experiment 1. This control experiment is part of Experiment 2 and should be treated alongside the other FPGA ambient EM and FPGA power-rail measurements. A companion PDF, `enf_based_chip_geolocation_additional_experiment_artifacts.pdf`, is included at the top level of the artifact package and contains the supplementary visual artifacts associated with this additional control experiment.

---

### 3.2 Experiment 3 – ENF Temporal Reliability and Replay Resilience

**Script:** `enf_analysis_top_reliability_stat.m`
**Data root:** `exp_inputs/TREND/`

**Goal:** Quantify how stable the ENF correlation is over **time-of-day, day-of-week, and week-to-week** at a fixed grid location. The same script also evaluates stale replay attacks by comparing selected target traces against older reference windows.

**Overview:**

1. **Inputs**
   - Walks the `TREND/WKxx/…` hierarchy described above.
   - For each (Week, Day, Time-of-Day, Trial) condition, loads:
     - `mains_pow_trace_ac.wav` – Ground-truth mains reference.
     - `fpga_em_trace_dc.wav` – Ambient EM trace near the board.

2. **ENF Extraction & Correlation**
   - Uses the same STFT+weighted-harmonic estimator as the pair-wise script.
   - Computes correlation **r** per file pair and applies the **Fisher transform** `z = atanh(r)` for parametric inference.

3. **Statistics**
   - Computes descriptive statistics for **r** and **z** (N, mean, SD, SEM, 95% CI).
   - Performs one-factor-at-a-time analyses on Fisher-z:
     - Time-of-Day (EMRN, MORN, AFTN, EVEN).
     - Day-of-Week (WED, THUR, SAT, SUN).
     - Week-to-Week (WK01 vs WK02) via paired t-tests on matched conditions.

4. **Visualization**
   - Dot-and-whisker (mean ± 95% CI) plots for the factors above.
   - Optional histograms and z-distribution diagnostics.

5. **Replay Attack Analysis**
   - Uses target traces configured in `config.replay_target_specs`.
   - Compares each selected target ambient trace against strictly older mains-reference windows.
   - Reports pooled authentic-vs-replay correlation statistics.
   - Exports `replay_attack_target_summary.csv` and a representative ENF overlay plot under `exp_results/replay_attack_resilience/`.

---

### 3.3 Experiment 4 – ENF Server-Room Robustness

**Script:** `enf_analysis_top_reliability_stat_server.m`
**Data root:** `exp_inputs/SRV_L/`

**Goal:** Evaluate how robust ENF extraction is inside a **high-density server-room environment**, focusing on:

- Sensor position on the server PSU (SBOX vs SPSU).
- Day-of-week effects (MON vs WED).

**Overview:**

1. **Inputs**
   - Traverses `SRV_L/SBOX/` and `SRV_L/SPSU/`, each with `MON/` and `WED/` subfolders.
   - Within each (Site, Day) folder, it loads:
     - `fpga_em_trace_dc.wav` – Ambient EM at the server PSU.
     - `mains_pow_trace_ac.wav` – Mains reference.
   - Five trials (`T01`–`T05`) per condition.

2. **ENF Extraction & Correlation**
   - STFT around 60 Hz harmonics.
   - Weighted-harmonic ENF estimator.
   - Pearson correlation **r** per trial, then Fisher-z for inference.

3. **Statistics**
   - Descriptive statistics on **r** and **z**.
   - One-way repeated-measures ANOVA on **z** for:
     - **Day**: MON vs WED (subjects = repeated trials, averaged over Site).
     - **Site**: SBOX vs SPSU (subjects = repeated trials, averaged over Day).
   - Reports effect sizes and 95% CIs, then back-transforms model estimates to **r**.

4. **Visualization**
   - Dot-and-whisker plots for Site and Day.
   - Boxplots and histograms for **r** and **z**.
   - Diagnostic z-distribution plots.

---

### 3.4 Experiment 5 – ENF Multi-Location Validation and Relay Attack Analysis

**Script:** `enf_analysis_top_multi_loc_stat.m`
**Data root:** `exp_inputs/MULTI/`

**Goal:** Quantify **cross-location separability** of ENF signatures by correlating a single local ambient-EM ENF against mains references drawn from multiple grid regions (e.g., US East / Texas / West 60 Hz, German 50 Hz). The script also evaluates relay-delay sensitivity by injecting controlled timing offsets into selected same-grid trace pairs.

**Overview:**

1. **Inputs**
   - Root directory contains trial folders (`T01/`, `T02/`, `T03/`, …) for each date.
   - Each `Txx` folder holds:
     - One local EM trace: `fpga_em_trace_dc_*lab.wav`.
     - Multiple mains references whose filenames encode **grid/region**:
       - `*_egrid_citya_lab.wav`, `*_egrid_citya_home.wav`, `*_egrid_worcester.wav`, `*_tgrid_richardson.wav`, `*_wgrid_tucson.wav`, `*_dresden.wav`, etc.
   - Subtrees `US_60/` and `DE_50/` correspond to 60 Hz and 50 Hz grids.

2. **ENF Extraction & Correlation**
   - STFT centered on the appropriate nominal frequency (50 or 60 Hz) and harmonics.
   - ENF estimation via weighted-harmonic approach for each trace.
   - For every EM trace, computes correlation **r** with every available mains reference in that trial folder.
   - Aggregates pairwise correlations into:
     - Per-folder lower-triangle correlation matrices.
     - Per-grid average correlation matrices and summaries.

3. **Statistics**
   - Descriptive statistics on **r** and **z = atanh(r)** for:
     - Intra-grid (same grid) vs inter-grid (cross-grid) comparisons.
   - Two-sample t-tests on **z** (intra-grid > inter-grid, Welch unequal-variance).
   - Optional parametric **Equal-Error Rate (EER)** estimates under a Gaussian model for the intra-grid vs inter-grid distributions, including “conservative” variants using reliability-derived variance.

4. **Relay Attack Analysis**
   - Set `config.ANALYSIS_MODE = 'relay_attack'` to run the relay-delay analysis.
   - Uses selected US_60 EGrid trace pairs.
   - Sweeps injected relay delays from −10 s to +10 s in 100 ms steps.
   - Reports the zero-delay baseline as the 0 injected-delay case.
   - Exports correlation-vs-delay plots such as `relay_attack_corr_vs_delay_us60_egrid_enf.svg`.

5. **Visualization**
   - Lower-triangle heatmaps with annotated correlation coefficients.
   - Per-grid average heatmaps (e.g., `US_60` and `DE_50` summaries).
   - Histograms and boxplots for **r** and **z**.
   - Optional plots of z-distributions and ROC/EER operating points.
   - Relay-delay correlation curves when running relay attack mode.

---

### 3.5 Pair-Wise Synthetic Forgery Analysis

**Script:** `enf_analysis_top_pair_wise_af.m`
**Data root:** one selected trace under `exp_inputs/MULTI/`

**Goal:** Evaluate whether synthetic ENF forgery and anti-forensics transformations can evade ENF-based checks. The script starts from an authentic AC mains trace, generates forged variants, and runs correlation and inter-harmonic consistency checks.

**Overview:**

1. **Inputs**
   - Loads one selected authentic trace, currently configured from `exp_inputs/MULTI/US_60/OCT/MON/T02/`.
   - The selected file is `mains_pow_trace_ac_egrid_citya_lab.wav` by default.

2. **Attack Variants**
   - **A0_2: ENF removal narrow cascade**
     - Applies a narrow-band cascade FIR bandstop operation around the ENF fundamental.
   - **A1_2: ENF noise fill narrow cascade**
     - Adds narrow-band noise fill matched to the authentic trace background spectrum after A0_2 removal.
   - **A2_2: synthetic ENF embedding narrow cascade**
     - Injects a peak-magnitude-matched FM synthetic ENF component into the A1_2 base.

3. **Detection**
   - **D0 ENF correlation check**
     - Extracts ENF from the authentic and forged traces, then reports authentic-vs-forged Pearson correlation.
   - **D1 inter-harmonic consistency check**
     - Compares ENF extracted at the fundamental and higher harmonics to detect inconsistency introduced by forgery.

4. **Outputs**
   - FFT harmonic attenuation report.
   - Per-attack FFT overlay plots.
   - Spectrogram comparison plots.
   - D0 ENF overlay plots.
   - D1 harmonic consistency plots.
   - Output directory: `exp_results/af_analysis_pairwise/`.

---

### 3.6 Core ENF Extraction Function

**File:** `proc_enf_analysis.p`

This compiled MATLAB function encapsulates the shared pipeline used by all top-level scripts:

1. Pre-processing & optional resampling.
2. STFT/spectrogram computation around the nominal ENF harmonics.
3. Weighted-harmonic ENF estimation per trace.
4. Temporal alignment and correlation computation.
5. Optional plot generation.

The top-level scripts configure and call this function; you do not need to invoke it directly.


## 4. Reproducing the Analysis

### 4.1 Requirements

To run the analysis end-to-end, you will need:

- **MATLAB R2018b or newer**.
- The **Signal Processing Toolbox**.
- Sufficient disk space to load the WAV files and write figures / logs.

No additional toolboxes are required beyond what the top-level scripts are already using.

### 4.2 Quick Start

1. **Clone or download this repository** to your local machine.

2. **Open MATLAB** and set the current folder to the repository root:

   ```matlab
   cd('PATH/TO/tches2027_enfchipgeoloc');
   ```

3. **Run scripts from the `exp_scripts/` directory**:

   ```matlab
   cd('exp_scripts');
   ```

4. **Run a pair-wise ENF feasibility analysis**:

   ```matlab
   run('enf_analysis_top_pair_wise.m');
   ```

5. **Run temporal reliability and replay-resilience analysis**:

   ```matlab
   run('enf_analysis_top_reliability_stat.m');
   ```

6. **Run server-room robustness analysis**:

   ```matlab
   run('enf_analysis_top_reliability_stat_server.m');
   ```

7. **Run multi-location baseline statistics**:

   ```matlab
   % In enf_analysis_top_multi_loc_stat.m:
   % config.ANALYSIS_MODE = 'baseline_stats';
   run('enf_analysis_top_multi_loc_stat.m');
   ```

8. **Run relay attack analysis**:

   ```matlab
   % In enf_analysis_top_multi_loc_stat.m:
   % config.ANALYSIS_MODE = 'relay_attack';
   run('enf_analysis_top_multi_loc_stat.m');
   ```

9. **Run pair-wise synthetic forgery / anti-forensics analysis**:

   ```matlab
   run('enf_analysis_top_pair_wise_af.m');
   ```

### 4.3 Notes on Selecting Pair-Wise Inputs

The pair-wise script analyzes one selected trace pair per run. Update the input variables near the top of `enf_analysis_top_pair_wise.m`:

- `file_1_path` – folder containing the AC mains reference.
- `file_2_path` – folder containing the sensed trace.
- `file_1_name` – mains reference file name without `.wav`.
- `file_2_name` – sensed trace file name without `.wav`.

For most datasets, `file_2_path = file_1_path`. For the FPGA datasets:

- Use `exp_inputs/FPGA/EM_TRACES/<BOARD>/DC_x/` with `fpga_em_trace_dc` for ambient EM analysis.
- Use `exp_inputs/FPGA/EM_TRACES/SAKU/ON_OFF_ON/` with `mains_pow_trace_ac` and `fpga_em_trace_dc` for the Experiment 2 ON-OFF-ON control trace.
- Use `exp_inputs/FPGA/POW_TRACES/<BOARD>/DC_x/` with `fpga_pow_trace_dc_a` or `fpga_pow_trace_dc_b` for board power-trace analysis.

For the ON-OFF-ON control trace, run `enf_analysis_top_pair_wise.m` with the `SAKU/ON_OFF_ON` path selection shown in the script comments. If chunk processing is enabled, the script can split the 200 s recording into the three control intervals (0-50 s, 50-120 s, and 120-200 s) for interval-wise correlation analysis.

For 50 Hz `MULTI/DE_50/` traces, update the nominal frequency settings in the pair-wise script from 60 Hz to 50 Hz before running.



