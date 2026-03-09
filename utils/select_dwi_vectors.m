function [adc_vec, d_vec, f_vec, dstar_vec] = select_dwi_vectors(data_vectors_gtvp, j, k, rpi, dwi_type)
% SELECT_DWI_VECTORS — Extract parameter vectors for the specified DWI pipeline.
%
%   Selects the appropriate ADC, D, f, and D* voxel vectors from the
%   data_vectors_gtvp struct array based on the DWI processing type.
%
%   Parameters
%   ----------
%   data_vectors_gtvp : struct array
%       3D struct array (patient x timepoint x repeat) containing voxel vectors.
%   j : integer
%       Patient index.
%   k : integer
%       Timepoint index.
%   rpi : integer
%       Repeat index.
%   dwi_type : integer (1, 2, or 3)
%       DWI processing type:
%         1 = Standard: raw DWI with conventional mono/bi-exponential fits
%         2 = dnCNN: DnCNN-denoised DWI signal with conventional fits
%         3 = IVIMnet: raw DWI with neural-network IVIM parameter estimation
%         Note: IVIMnet reuses standard ADC (mono-exponential fit is model-free
%         and unaffected by the IVIM fitting method choice).
%
%   Returns
%   -------
%   adc_vec : numeric vector
%       ADC values for each voxel.
%   d_vec : numeric vector
%       Diffusion coefficient (D) values.
%   f_vec : numeric vector
%       Perfusion fraction (f) values.
%   dstar_vec : numeric vector
%       Pseudo-diffusion coefficient (D*) values.
%
%   See also: compute_summary_metrics, compute_adc_metrics, compute_ivim_metrics

    switch dwi_type
        case 1
            adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector;
            d_vec     = data_vectors_gtvp(j,k,rpi).d_vector;
            f_vec     = data_vectors_gtvp(j,k,rpi).f_vector;
            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector;
        case 2
            adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector_dncnn;
            d_vec     = data_vectors_gtvp(j,k,rpi).d_vector_dncnn;
            f_vec     = data_vectors_gtvp(j,k,rpi).f_vector_dncnn;
            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_dncnn;
        case 3
            adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector;
            d_vec     = data_vectors_gtvp(j,k,rpi).d_vector_ivimnet;
            f_vec     = data_vectors_gtvp(j,k,rpi).f_vector_ivimnet;
            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_ivimnet;
    end
end
