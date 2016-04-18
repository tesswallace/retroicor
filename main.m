% Script to perform modified RETROICOR correction on BOLD image data
% Compiled by Tess Wallace (Department of Radiology, University of Cambridge)

% -------------------------------------------------------------------------
% Load DICOM images
% -------------------------------------------------------------------------
dicom_dir = 'E:\DICOM\Behold_Vol_030_20150327\Mri_Breast_Bilateral_W_Contrast - 0\BOLD_R_O2_Carb_10\';
[s_bold_t_series, dicom_header] = read_dicom_dir(dicom_dir);

% -------------------------------------------------------------------------
% Draw ROI for Analysis
% -------------------------------------------------------------------------
[ BW ] = draw_roi(mean(s_bold_t_series,3));

% -------------------------------------------------------------------------
% Physiological Data Directories
% -------------------------------------------------------------------------
PPGfilename = 'E:\Data\TraceData\Behold_Vol_030\PPG_Resp_Data\PPGData_Oxygen_Carbogen';
RESPfilename = 'E:\Data\TraceData\Behold_Vol_030\PPG_Resp_Data\RESPData_Oxygen_Carbogen';

% -------------------------------------------------------------------------
% Perform RETROICOR
% -------------------------------------------------------------------------
% Define start of scan relative to physiological data acquisition
% GE: calculated as (length scan - 200 ms + TE)
start_time = 34630; % ms 

% Fifth order Fourier series calculated to model cardiac (C) and
% respiratory (R) phases
% Four multiplicative (X) terms calculated
[design_matrix, TR_phs] = mod_retroicor(PPGfilename, RESPfilename, start_time, dicom_header);

% -------------------------------------------------------------------------
% Perform Linear Regression using GLM
% -------------------------------------------------------------------------
HW_flag = 1; % 1: linear, 2: quadratic polynomial to correct for drift
retroicor_flag = 1; % 0: off, 1: on

% Specify order of RETROICOR correction e.g. 2R2C1X
order.C = 1; % nth order cardiac (C=0-5)
order.R = 2; % nth order respiratory (R=0-5)
order.X = 1; % n multiplicative terms (X=0-4)

[s_bold_t_series_phys, yDelta, adjRsq] = glm_regression(s_bold_t_series, design_matrix, order, dicom_header, BW, HW_flag, retroicor_flag);