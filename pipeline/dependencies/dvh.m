function [dvhparams, dvh_values] = dvh(dose_map, struct, dims, nbins, dvhparams_in, normalize)
%dvh A function to calculate the cumulative dose-volume histogram (DVH) and
%   DVH parameters from a dose map and structure mask.
%
%   The primary purpose of this function is to calculate DVH parameters
%   (e.g. D99%, V40Gy, D0.5cc, ...) given a dose map and binary structure
%   mask.
%   If desired, the secondary output of this function is the cumulative DVH
%   as a nbins x 2 matrix, where the first column is the dose values of the
%   bins (i.e. the DVH x-axis) and the second column is the (by default
%   normalized) volume of the structure receiving at least that amount of
%   dose (i.e. DVH y-axis). The DVH can then be simply plotted with
%   plot(dvh_values(:, 1), dvh_values(:, 2))
% 
%   Inputs:
%   - dose_map: 2D/3D array with the dose values
%   - struct: logical array of the structure mask
%   - dims: 1 x ndims vector with a voxel sizes in cm(!) in each dimension
%   - nbins (optional, default 2000): number of DVH bins, maximum bin value
%       is always the maximum value in dose_map
%   - dvhparams_in (optional): name-value arguments specifying the DVH
%       parameters to be calculated. Syntax:
%      - 'Dperc' (%): reports the minimal dose (Gy) received by the given volume percentage
%      - 'Dvol' (cc): reports the minimal dose (Gy) received by the given volume
%      - 'Vperc' (Gy): reports the volume (%) receiving at least the given dose
%      - 'Vvol' (Gy): reports the volume (cc) receiving at least the given dose
%   - 'Normalize' (optional, default true): return the DVH volumes scaled
%       between 0 and 100% (true), or in cc (false)
% 
%   Outputs:
%   - dvhparams: table with the calculated DVH parameters
%   - dvh_values: nbins x 2 matrix, with first column the dose bins and
%       second column the normalized or absolute volume.
%   
%   Examples:
%   >> params = dvh(dose_map, gtv_mask, [0.3 0.3 0.3], 'Dperc', 99, 'Vvol', 45)
%    params =
% 
%      1×2 table
% 
%        D99% (Gy)    V45.0cc (cc)
%        _________    ____________
% 
%         40.711         29.889   

%
%
%   20230414 G. Grimbergen, Department of Radiation Oncology, UMC Utrecht

arguments 
    dose_map
    struct logical
    dims (1, :)
    nbins double = 2000
    dvhparams_in.Dperc
    dvhparams_in.Dvol
    dvhparams_in.Vperc
    dvhparams_in.Vvol
    normalize.Normalize = true
end

% Make a vector with the dose bins
dose_bins = linspace(0, max(dose_map, [], 'all'), nbins);

% Get a vector of all voxels within the structure
voxels = dose_map(struct);

% The following part is vectorized for better performance. This is way
% faster than looping over all bins and doing a > operation for every bin
% value
% Create a voxels x bins matrix
voxels_mat = repmat(voxels, 1, length(dose_bins));

% Compare every column against the dose bins to get a binary matrix
diff_mat = voxels_mat > dose_bins;

% The number of 1's in a column is the amount of voxels that is larger
% than a certain dose bin
volume = sum(diff_mat);

volume_normalized = 100 * (volume ./ sum(struct, "all")); % In percentage
volume = volume * prod(dims); % Convert number of voxels to volume

% There is probably a better way to parse a binary input like this than
% another name-value pair, but I'm still finding my way around Matlab's
% argument block functionality so bear with me here.
if normalize.Normalize
    dvh_values = [dose_bins; volume_normalized]';
else
    dvh_values = [dose_bins; volume]';
end


% If there are any, calculate the requested DVH parameters and put them in
% a nice Matlab-friendly table (that also allows %-signs in the variable
% names, as opposed to structure arrays).
dvhparams = table();
params = fieldnames(dvhparams_in);
for i = 1:length(params)
    current_param = params{i};
    switch current_param
        case 'Dperc'
            [~, vol_idx]=min(abs(volume_normalized - dvhparams_in.(current_param)));
            parameter_name = sprintf('D%i%% (Gy)', dvhparams_in.(current_param));
            dvhparams.(parameter_name) = dose_bins(vol_idx);
        case 'Dvol'
            [~, vol_idx]=min(abs(volume - dvhparams_in.(current_param)));
            parameter_name = sprintf('D%2.1fcc (Gy)', dvhparams_in.(current_param));
            dvhparams.(parameter_name) = dose_bins(vol_idx);
        case 'Vperc'
            [~, dose_idx]=min(abs(dose_bins - dvhparams_in.(current_param)));
            parameter_name = sprintf('V%iGy (%%)', dvhparams_in.(current_param));
            dvhparams.(parameter_name) = volume_normalized(dose_idx);
        case 'Vvol'
            [~, dose_idx]=min(abs(dose_bins - dvhparams_in.(current_param)));
            parameter_name = sprintf('V%2.1fGy (cc)', dvhparams_in.(current_param));
            dvhparams.(parameter_name) = volume(dose_idx);
    end
end



end