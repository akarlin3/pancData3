function adc = fit_adc_mono(dwi,bvals)
% simple monoexponential fit of 3D dwi scan
% format:   dwi [Ny,Nx,Nz,nbvals)
%           bvals: [nbvals,1] (first one should be zero)

dwi_x = size(dwi,2);
dwi_y = size(dwi,1);
dwi_z = size(dwi,3);
tmp0 = reshape(dwi,[dwi_y*dwi_x*dwi_z,length(bvals)]);
adc_tmp = -bvals(2:end) \ permute(log(tmp0(:,2:end)./tmp0(:,1)),[2 1]);
adc = abs(reshape(adc_tmp,[dwi_y, dwi_x, dwi_z]));