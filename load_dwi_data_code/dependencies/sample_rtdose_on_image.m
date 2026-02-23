function dose_sampled = sample_rtdose_on_image(imagelocation,rtdosefile)
% inputs: imagelocation: either directory containing dicom files, or list of files
%                rtdose: rtdose dicom file * note, this assumes isoropic dose
%                        sampling (i.e. x res = y res = z res)

%%%%%%% load three-dimensional image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[image,spatial,dim] = dicomreadVolume(imagelocation);
image = squeeze(image);
image = image - 3614;
image_origin = spatial.PatientPositions(1,:);
image_spacing = spatial.PixelSpacings(1,:);
image_spacing(3) = spatial.PatientPositions(2,3) - spatial.PatientPositions(1,3);
image_size = size(image);
% coordinates of image
x_image = zeros(image_size(1),1);
y_image = zeros(image_size(2),1);
z_image = zeros(image_size(3),1);
for kk = 1:image_size(1)
    x_image(kk) = image_origin(1) + image_spacing(1)*(kk-1);
end
for ll = 1:image_size(2)
    y_image(ll) = image_origin(2) + image_spacing(2)*(ll-1);
end
for mm = 1:image_size(3)
    z_image(mm) = image_origin(3) + image_spacing(3)*(mm-1);
end


%%%%% read DICOM RT dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtdose_info = dicominfo(rtdosefile);
rtdose_data = dicomread(rtdose_info);
rtdose_data = squeeze(rtdose_data);
rtdose_origin = rtdose_info.ImagePositionPatient;
rtdose_spacing(1:2) = rtdose_info.PixelSpacing;

[xxx_image,yyy_image,zzz_image] = meshgrid(x_image,y_image,z_image);

%%
rtdose_spacing(3) = rtdose_spacing(1); %rtdose_info.SliceThickness;
rtdose_size = size(rtdose_data);
% to convert raw data to dose
rtdose_gridscaling = rtdose_info.DoseGridScaling;
rtdose = rtdose_gridscaling*double(rtdose_data);
% coordinates of image
x_rtdose = zeros(rtdose_size(2),1);
y_rtdose = zeros(rtdose_size(1),1);
z_rtdose = zeros(rtdose_size(3),1);
for kk = 1:rtdose_size(2)
    x_rtdose(kk) = rtdose_origin(1) + rtdose_spacing(1)*(kk-1);
end
for ll = 1:rtdose_size(1)
    y_rtdose(ll) = rtdose_origin(2) + rtdose_spacing(2)*(ll-1);
end
for mm = 1:rtdose_size(3)
    z_rtdose(mm) = rtdose_origin(3) + rtdose_spacing(3)*(mm-1);
end
[xxx_rtdose,yyy_rtdose,zzz_rtdose] = meshgrid(x_rtdose,y_rtdose,z_rtdose);

%%

dose_sampled = interp3(xxx_rtdose, yyy_rtdose, zzz_rtdose, rtdose, xxx_image, yyy_image, zzz_image);

end