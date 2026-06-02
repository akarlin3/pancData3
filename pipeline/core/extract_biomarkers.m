function bio = extract_biomarkers(mask, maps, meta, dncnn_maps, dncnn_mask, ivimnet_maps)
    % EXTRACT_BIOMARKERS Extract voxel-level parameter vectors within a mask
    %   Extracted local helper from process_single_scan.m. Returns a scan
    %   struct (init_scan_structs template) populated with ADC/D/f/D* vectors
    %   from MAPS inside MASK, scan metadata from META, and optional
    %   DnCNN/IVIMnet vectors when those maps are present.
    [tmp, ~] = init_scan_structs(1, 1);
    bio = tmp;
    mask_idx = (mask == 1);
    bio.adc_vector = maps.adc_map(mask_idx);
    bio.d_vector   = maps.d_map(mask_idx);
    bio.f_vector   = maps.f_map(mask_idx);
    bio.dstar_vector = maps.dstar_map(mask_idx);
    bio.ID = meta.id;
    bio.MRN = meta.mrn;
    bio.LF = meta.lf;
    bio.Immuno = meta.immuno;
    bio.Fraction = meta.fi;
    bio.Repeatability_index = meta.rpi;
    bio.vox_vol = meta.vox_vol;
    bio.vox_dims = meta.vox_dims;
    if isfield(dncnn_maps, 'd_map_dncnn') && ~isempty(dncnn_maps.d_map_dncnn)
        dm = (dncnn_mask == 1);
        bio.d_vector_dncnn     = dncnn_maps.d_map_dncnn(dm);
        bio.f_vector_dncnn     = dncnn_maps.f_map_dncnn(dm);
        bio.dstar_vector_dncnn = dncnn_maps.dstar_map_dncnn(dm);
        bio.adc_vector_dncnn   = dncnn_maps.adc_map_dncnn(dm);
    end
    if isfield(ivimnet_maps, 'D_ivimnet') && ~isempty(ivimnet_maps.D_ivimnet)
        bio.d_vector_ivimnet     = ivimnet_maps.D_ivimnet(mask_idx);
        bio.f_vector_ivimnet     = ivimnet_maps.f_ivimnet(mask_idx);
        bio.dstar_vector_ivimnet = ivimnet_maps.Dstar_ivimnet(mask_idx);
    end
end
