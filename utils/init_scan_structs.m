function [gtvp_structs, gtvn_structs] = init_scan_structs(n_fx, n_rp)
% INIT_SCAN_STRUCTS Create empty GTVp/GTVn struct arrays with consistent field ordering.
%   [gtvp_structs, gtvn_structs] = init_scan_structs(n_fx, n_rp) returns two
%   n_fx-by-n_rp struct arrays pre-populated with all biomarker, dosimetry,
%   clinical, and deep-learning fields set to []. This ensures a uniform
%   field order across the pipeline and prevents struct concatenation errors
%   when fields are added in varying order.
%
%   Used by load_dwi_data and process_single_scan to initialise storage
%   before populating per-scan results.
%
%   Inputs:
%       n_fx  — number of fractions (rows)
%       n_rp  — number of repeatability indices (columns)
%
%   Outputs:
%       gtvp_structs — n_fx x n_rp struct array for primary GTV data
%       gtvn_structs — n_fx x n_rp struct array for nodal GTV data

    % --- Field Definitions ---
    % Each scan produces voxel-level vectors of diffusion and dosimetric
    % parameters extracted from within the GTV mask. The struct fields are
    % organized into logical groups:
    %
    % IVIM / ADC Biomarker Vectors (Standard fitting):
    %   adc_vector   - Apparent Diffusion Coefficient from mono-exponential fit
    %                  (units: mm^2/s). Reflects overall water mobility in tissue;
    %                  decreases with high cellularity (viable tumor).
    %   d_vector     - True diffusion coefficient from bi-exponential IVIM fit
    %                  (units: mm^2/s). Separates pure tissue diffusion from
    %                  perfusion effects, providing a purer cellularity measure.
    %   f_vector     - Perfusion fraction from IVIM fit (dimensionless, 0-1).
    %                  Represents the fraction of signal arising from blood flow
    %                  in the capillary network (pseudo-diffusion compartment).
    %   dstar_vector - Pseudo-diffusion coefficient from IVIM fit (units: mm^2/s).
    %                  Reflects microvascular blood flow velocity; typically 10x
    %                  larger than D and highly variable in pancreatic tissue.
    %
    % Dosimetry Fields:
    %   dose_vector  - Per-voxel radiation dose (Gy) within the GTV, resampled
    %                  from the RT dose grid onto the DWI coordinate space.
    %   dvh          - Dose-volume histogram for the GTV subvolume.
    %   d95          - Dose covering 95% of the GTV volume (Gy); a standard
    %                  radiotherapy quality metric.
    %   v50gy        - Fraction of GTV receiving >= 50 Gy; used to assess
    %                  dose coverage adequacy.
    %
    % Deep Learning Denoised Variants:
    %   *_dncnn      - IVIM/ADC parameters fitted after DnCNN denoising.
    %                  DnCNN reduces noise in low-SNR DWI images (especially at
    %                  high b-values), improving IVIM fit stability.
    %   *_ivimnet    - IVIM parameters from IVIMnet (neural network-based IVIM
    %                  fitting). Bypasses conventional least-squares fitting,
    %                  which is unstable for the bi-exponential IVIM model.
    %
    % Clinical / Metadata Fields:
    %   ID                  - Patient study identifier
    %   MRN                 - Medical Record Number (de-identified in outputs)
    %   LF                  - Local failure flag (0=LC, 1=LF, 2=competing risk)
    %   Immuno              - Immunotherapy flag (treatment stratification)
    %   Fraction            - Treatment fraction number (longitudinal index)
    %   Repeatability_index - Scan repeat index within the fraction
    %   vox_vol             - Voxel volume in mm^3 (from DICOM header); needed
    %                         to convert voxel counts to physical volumes for
    %                         DVH computation and tumor volume tracking.
    empty_entry = struct( ...
        'adc_vector', [], 'd_vector', [], 'f_vector', [], 'dstar_vector', [], ...
        'dose_vector', [], 'dvh', [], 'd95', [], 'v50gy', [], ...
        'd_vector_dncnn', [], 'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
        'adc_vector_dncnn', [], ...
        'd_vector_ivimnet', [], 'f_vector_ivimnet', [], 'dstar_vector_ivimnet', [], ...
        'ID', [], 'MRN', [], 'LF', [], 'Immuno', [], ...
        'Fraction', [], 'Repeatability_index', [], 'vox_vol', []);

    % --- Pre-allocate Struct Arrays ---
    % Using repmat with a fully-defined template struct guarantees that all
    % entries share identical field names in identical order. This is critical
    % because MATLAB struct array concatenation (e.g., [s1; s2]) fails if the
    % structs have different field orderings. Since process_single_scan
    % populates these structs inside a parfor loop (where execution order is
    % non-deterministic), consistent initialization prevents intermittent
    % concatenation failures that would be difficult to reproduce.
    %
    % The n_fx x n_rp layout mirrors the clinical data structure:
    %   Rows    = treatment fractions (e.g., baseline, mid-treatment, end)
    %   Columns = repeat scans within each fraction (for reproducibility analysis)
    gtvp_structs = repmat(empty_entry, n_fx, n_rp);
    gtvn_structs = repmat(empty_entry, n_fx, n_rp);
end
