function prev_states = suppress_core_warnings()
% SUPPRESS_CORE_WARNINGS — Suppress expected warnings from extract_tumor_core.
%
%   prev_states = suppress_core_warnings()
%
%   Returns a struct array of previous warning states so the caller can
%   restore them after the operation (e.g., via onCleanup).
%
%   These warnings are expected during batch core method comparison and
%   multi-method dosimetry, where some methods are intentionally run on
%   data that does not meet their prerequisites (e.g., spectral clustering
%   with too few voxels, active contours without 3D mask).

    ids = {
        'extract_tumor_core:tooFewForSpectral'
        'extract_tumor_core:no3DForActiveContours'
        'extract_tumor_core:no3DForRegionGrowing'
        'extract_tumor_core:fdmBaseline'
        'extract_tumor_core:fdmNoBaseline'
        'extract_tumor_core:noDValues'
        'extract_tumor_core:noIVIMValues'
        'extract_tumor_core:noSpectralCluster'
    };

    prev_states = struct('identifier', ids, 'state', '');
    for i = 1:numel(ids)
        prev_states(i) = warning('off', ids{i});
    end
end
