function [design_matrix, TR_phs] = mod_retroicor(PPGfilename, RESPfilename, start_time, dicom_header)

% Matlab implementation of RETROICOR correction algorithm
% Adapted from http://cbi.nyu.edu/software/

% Function inputs:
    % PPGfilename: pointer to file containing recorded PPG data
    % RESPfilename: pointer to file containing recorded RESP data
    % start_time: start time of scan relative to beginning of physiological 
    % data recording
    % dicom_header: struct containing relevant DICOM header information
    
% Function outputs:
    % design_matrix: matrix of cardiac (C), respiratory (R) and 
    % multiplicative (X) sine and cosine terms, based on the C and R phases
    % at each image time point
    % TR_phs: C and R phases at each image time point 

order = 5; % nth order correction

% -------------------------------------------------------------------------
% Read in physiological data
% -------------------------------------------------------------------------

fileID = fopen(PPGfilename, 'r');
PPGData = fscanf(fileID, '%f');
fclose(fileID);

fileID = fopen(RESPfilename, 'r');
RESPData = fscanf(fileID, '%f');
fclose(fileID);

% Define sample rate
PPGsr = 10; %sample rate in ms
RESPsr = 40; %sample rate in ms

tPPG = PPGsr*(1:1:length(PPGData));
tRESP = RESPsr*(1:1:length(RESPData));

h1=figure;
subplot(2,1,1);
plot(tPPG(1:2000), PPGData(1:2000));
title('PPG Trace');
xlabel('Time (ms)');
subplot(2,1,2);
plot(tRESP(1:500), RESPData(1:500));
title('RESP Trace');
xlabel('Time (ms)');

dt = 1/40;  

%Find TR times
TR_len = dicom_header.TR/PPGsr; % ms
numTR = dicom_header.tdim;
TR_start = start_time/PPGsr + (0:TR_len:((numTR-1)*TR_len)); %in samples

%% ------------------------------------------------------------------------
% Compute Cardiac Phase
% -------------------------------------------------------------------------

h2=figure;
subplot(2,1,1);
plot(1:1:length(PPGData), PPGData);
hold;
pulse_axis = axis;
for i=1:numTR
    plot([TR_start(i) TR_start(i)],pulse_axis(3:4),'g');
end
    
smp = 1:1:length(PPGData);

% Find peaks in pulse
lm = local_max(PPGData); 
c_peaks = smp(lm);
plot(c_peaks,PPGData(c_peaks),'r.');

const=input('Enter threshold for peak detection: ','s');
const = str2double(const);

pulse_mean = mean(PPGData);
plot([0 length(PPGData)],[const+pulse_mean const+pulse_mean],'y');
rmpk1 = find(PPGData(c_peaks) < pulse_mean+const); 
plot(c_peaks(rmpk1),PPGData(c_peaks(rmpk1)),'k.');
    
%manual peak editing
medit=input('Manual Edit (y/n)? ','s');

if medit == 'y'
    dch=datacursormode;
    set(dch,'DisplayStyle','datatip','SnapToDataVertex','on');
    mp_count=0; getpoints=1;
    while getpoints
        cont=1;
        cont=input('Click on a datapoint to remove, then hit enter. Type 0 when done:');
        if ~cont
            disp('Done!');
            getpoints = 0; 
        else
            mp_count=mp_count+1;
            dcinfo = getCursorInfo(dch);
            man_points(mp_count) = dcinfo.DataIndex;
            delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
            plot(c_peaks(man_points(mp_count)),PPGData(c_peaks(man_points(mp_count))),'c.');
        end 
    end
    delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
    plot(c_peaks(man_points),PPGData(c_peaks(man_points)),'k.');
    rmpk = union(rmpk1,man_points);
else
    rmpk = rmpk1;
end

c_peaks = setdiff(c_peaks,c_peaks(rmpk));

for i=1:length(PPGData)
    if i < c_peaks(1)
        c_phs(i) = NaN;
    elseif i >= c_peaks(end)
        c_phs(i) = NaN;
    else
        prev_peak = max(find(c_peaks <=i));
        t1 = c_peaks(prev_peak);
        t2 = c_peaks(prev_peak+1);
        c_phs(i) = 2*pi*(i - t1)/(t2-t1); %find cardiac phase for each acquisition
    end
end

plot(c_phs,'m');
xlabel('Samples');
title('Computing Cardiac Phase');
ylabel('mV');

TR_phs(:,1) = c_phs(TR_start);
subplot(2,1,2);
scatter(TR_phs(:,1),PPGData(TR_start));
xlabel('Phase in Cardiac Cycle (Radians)');
ylabel('mV');

% Fit nth order fourier series to estimate phase
for i = 1:order
   dm_c_phs(:,(i*2)-1) = cos(i*TR_phs(:,1));
   dm_c_phs(:,i*2) = sin(i*TR_phs(:,1));
end

h3=figure;
imagesc(dm_c_phs);
colormap(gray);
title('Design Matrix (C)');

%% ------------------------------------------------------------------------
% Compute Respiratory Phase
% -------------------------------------------------------------------------

% Normalize to range of 0 to 1
resp_range = range(RESPData);
resp_norm = (RESPData - min(RESPData) ) / resp_range;

h4=figure;
subplot(3,1,1);
plot(resp_norm);
hold;
 
% Histogram-equalized transfer function between respiratory amplitude and
%resp phase
nbins = 100;
[resp_hist,bins] = hist(resp_norm,nbins);
resp_transfer_func = [0 (cumsum(resp_hist) / sum(resp_hist))];
kern_size = 1/dt - 1; 
resp_smooth = conv(resp_norm,ones(kern_size,1),'same'); %smoothed version for taking derivative
resp_diff = [diff(resp_smooth);0]; %derivative dR/dt
r_phs = pi*resp_transfer_func(round(resp_norm * nbins)+1)' .* sign(resp_diff);

plot(resp_smooth / max(resp_smooth),'g'); %plot smoothed version
axis([0 length(resp_norm) 0 1]);

subplot(3,1,2);
plot(r_phs,'m');
hold;
 for i=1:numTR
    plot([TR_start(i)/4 TR_start(i)/4],[-pi pi],'g');
 end
axis([0 length(r_phs) -pi pi]);

% Get TR phase
TR_phs(:,2) = r_phs(round(TR_start/4));
subplot(3,1,3);
scatter(TR_phs(:,2),resp_norm(round(TR_start/4)));
axis([-pi pi 0 1]);
xlabel('Phase in Respiratory Cycle (Radians)');
ylabel('Normalized Respiration Belt (norm V) ');

% Fit nth order fourier series to estimate phase
for i = 1:order
    dm_r_phs(:,(i*2)-1) = cos(i*TR_phs(:,2));
    dm_r_phs(:,i*2) = sin(i*TR_phs(:,2));
end
    
h5=figure;
imagesc(dm_r_phs);
colormap(gray);
title('Design Matrix (R)');

%% ------------------------------------------------------------------------
% Compute Multiplicative Terms
% -------------------------------------------------------------------------

i=1;
for C = 1:2
    for D = 1:2
    dm_cr_phs(:,(i*4)-3) = sin(C*TR_phs(:,1)+D*TR_phs(:,2));
    dm_cr_phs(:,(i*4)-2) = cos(C*TR_phs(:,1)+D*TR_phs(:,2));
    dm_cr_phs(:,(i*4)-1) = sin(C*TR_phs(:,1)-D*TR_phs(:,2));
    dm_cr_phs(:,(i*4)) = cos(C*TR_phs(:,1)-D*TR_phs(:,2));
    i=i+1;
    end
end

h6=figure;
imagesc(dm_cr_phs);
colormap(gray);
title('Design Matrix (X)');

%% ------------------------------------------------------------------------
% Assign output design matrix
% -------------------------------------------------------------------------

design_matrix = [dm_c_phs dm_r_phs dm_cr_phs];

fin=input('Hit Enter to close plots and exit');
close([h1 h2 h3 h4 h5 h6]);

