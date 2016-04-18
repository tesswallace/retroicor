function [s_bold_t_series, dicom_header] = read_dicom_dir(dicom_dir)

% This function reads in a dicom series from a specified directory
% and stores the image volume in a variable s_bold_t_series, and key
% information in a struct dicom_header

filter = '*.dcm';

dcmseries = dir(strcat(dicom_dir, filter));

filename = strcat(dicom_dir, dcmseries(1).name);

info = dicominfo(filename);

dicom_header.xdim=info.Width;
dicom_header.ydim=info.Height;
dicom_header.zdim = 1;
dicom_header.tdim = info.NumberOfTemporalPositions;
dicom_header.PixelSpacing = info.PixelSpacing;
dicom_header.SliceThickness = info.SliceThickness;
dicom_header.TR = info.RepetitionTime;

s_bold_t_series = zeros(dicom_header.ydim, dicom_header.xdim, dicom_header.tdim);

for i=1:dicom_header.tdim
    
    filename = strcat(dicom_dir, dcmseries(i).name);
    
    s_bold_t_series(:,:,i) = dicomread(filename);

end

s_bold_t_series = double(s_bold_t_series);

end

