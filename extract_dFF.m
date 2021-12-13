function [dFF465_M, dFF565_M] = extract_dFF(raw_data, Fs, exp_ID, start_time, figtoggle)
%% DeltaF_highpass Smoothes and Processes Data for Fiber Photometry
% Converts raw data to dF/F values usinges a series of filters to 
% remove artifacts, de-noise traces, and correct for photobleaching and
% movement. The baseline is estimated using a lowpass filter with a low
% frequency cutoff and dF/F is returned as a percentage. Most inputs to
% this function can be extracted using "extract_TDT_Tank.m"
% James Maksymetz and Max Joffe April 2020.
% 
% https://github.com/ThomasAkam/photometry_preprocessing

%% Inputs
%
% # raw_data, array containing continuous data (seconds and mV)
% (row 1:time, 2: 465, 3:405, 4:565)
% # Fs, final sampling frequency (accounting for downsampling)
% # exp_ID, name to be appended to files and figures
% # start_time, start of experiment, (same units as raw_data; seconds)
% # Toggle for figures (1 = make)

%% Output
%
% # dFF465_M is a 3-column matrix containing times, dF/F 465, dF/F 405
% # dFF565_M is a 3-column matrix containing times, dF/F 565, dF/F 405
%
% dFF565_M returns an empty matrix if Data array only has 3 rows

%% Polymorphic Options
if nargin == 1
    error('Not enough input variables to proceed.')
else
    if nargin < 5
        figtoggle = 1;
        
        if nargin < 4
            start_time = 0;
        
        if nargin < 3
            exp_ID = char(datetime('now')); % timestamp will replace Block if not provided
            exp_ID = exp_ID(exp_ID~=':'); % remove colons
        
            if nargin == 0
                figtoggle = 1; 
                BlockDir=uigetdir(cd, 'Select TDT Photometry Block'); % get folder name
                cd(BlockDir); % change directory to photometry folder
                addpath(BlockDir);
                [raw_data, Fs, exp_ID] = extract_TDT_Tank(BlockDir); % extract raw data from Tank
            end
        end                
        end
    end
end



%% Signal Pre-Processing

% Trim beginning of file to avoid fitting pre-experiment traces
if start_time ~= 0
    raw_data = raw_data(:,floor(start_time*Fs):end);
end
    
% Separate data table into at least 3 column vectors
Ts = double(raw_data(1,:)');
Ch465 = double(raw_data(2,:)');
Ch405 = double(raw_data(3,:)');

[num_sig, ~] = size(raw_data);

if num_sig == 4
    Ch565 = double(raw_data(4,:)');
end

%% Median Filtering to Remove Artifacts
dn465 = medfilt1(Ch465, 1, [], 5);
dn405 = medfilt1(Ch405, 1, [], 5);

%% Low Pass Filter to Remove Noise
fc = 10; % lowpass filter value, (default = 10Hz)
[b,a] = butter(2, fc/(Fs/2), 'low');
dn465 = filtfilt(b, a, dn465);
dn405 = filtfilt(b, a, dn405);

if num_sig == 4
    dn565 = medfilt1(Ch565, 1, [], 5);
    dn565 = filtfilt(b, a, dn565);
end
 
% Plot raw signals pre- and post-filtering
if figtoggle == 1
f3 = figure;

if num_sig == 3
subplot (2,1,1); 
elseif num_sig == 4
    subplot (3,1,1);
end    
plot(Ts, Ch465, 'k');
hold on;
plot(Ts, dn465, 'g');
title('Raw and Denoised/Filtered Signals');
ylabel('mV (465 nm)');
legend('Raw Signal', sprintf('Med/Low Filt Signal, Cutoff %d Hz', fc));
axis tight

if num_sig == 4
    subplot (3,1,2);    
plot(Ts, Ch565, 'k');
hold on;
plot(Ts, dn565, 'r');
ylabel('mV (565 nm)');
legend('Raw Signal', sprintf('Med/Low Filt Signal, Cutoff %d Hz', fc));
axis tight
f3.Position = [500 100 550 700];
end

if num_sig == 3
subplot (2,1,2);
elseif num_sig == 4
    subplot(3,1,3);
end
plot(Ts, Ch405, 'k');
hold on;
plot(Ts, dn405, 'm');
ylabel('mV (405 nm)');
xlabel('Time (s)');
legend('Raw Signal', sprintf('Med/Low Filt Signal, Cutoff %d Hz', fc),...
       'Location','Best');
axis tight

savefig(strcat(exp_ID,' raw and filtered'));
saveas(f3,strcat(exp_ID,' raw and filtered.jpg'));
end

%% Low-cutoff Highpass Filter to Correct for Photobleaching
fc = 0.001; % cutoff of the highpass filter (default = 0.001 Hz, 16 min)
[b,a] = butter(2, fc/(Fs/2), 'high');
hf465 = filtfilt(b, a, dn465);
hf405 = filtfilt(b, a, dn405);

if num_sig == 4
    hf565 = filtfilt(b, a, dn565);
end

% Plot signals corrected for photobleaching
%if figtoggle == 1
%f4=figure;
%plot(Ts, hf465, 'g');
%hold on;
%plot(Ts, hf405-20, 'm'); % 405 signal shifted down for visualuzation purposes
%title('Photobleaching Corrected Signals by Highpass Filter');
%ylabel('mV');
%xlabel('Time (s)');
%if num_sig == 4
%    plot(Ts, hf565, 'r');
%    legend('465 nm', '565 nm', '405 nm (down 20 mV)');
%elseif num_sig == 3
%    legend('465 nm', '405 nm (down 20 mV)', 'Location','Best');
%end
%axis tight
%savefig(strcat(exp_ID,' photobleaching corrected'));
%saveas(f4,strcat(exp_ID,' photobleaching corrected.jpg'));
%end

%% Estimate Baseline F with a Low-Cutoff Lowpass Filter
fc = 0.001;
[b,a] = butter(2, fc/(Fs/2), 'low');
baseline465 = filtfilt(b, a, dn465);
baseline405 = filtfilt(b, a, dn405);

if num_sig == 4
    baseline565 = filtfilt(b, a, dn565);
end

% Plot the baseline fluorescence
if figtoggle == 1
f5 = figure;

if num_sig == 3
subplot(2,1,1);
elseif num_sig == 4
    subplot (3,1,1);
end
plot(Ts, hf465, 'g');
hold on;
plot(Ts, baseline465 - mean(baseline465), 'k');
title('Baseline Fluorescence');
ylabel('mV');
legend('Filtered 465nm', ['Baseline - ', num2str(mean(baseline465))], 'Location','Best');
axis tight

if num_sig == 4
    subplot (3, 1, 2)
plot(Ts, hf565, 'r');
hold on;
plot(Ts, baseline565 - mean(baseline565), 'k');
ylabel('mV');
legend('Filtered 565nm', ['Baseline - ', num2str(mean(baseline565))]);
axis tight
f5.Position = [500 100 550 700];
end

if num_sig == 3
subplot(2,1,2);
elseif num_sig == 4
    subplot(3,1,3);
end
plot(Ts, hf405, 'm');
hold on;
plot(Ts, baseline405 - mean(baseline405), 'k');
ylabel('mV');
xlabel('Time (s)');
legend('Filtered 405nm', ['Baseline - ', num2str(mean(baseline405))], 'Location','Best');
axis tight

savefig(strcat(exp_ID,' processed'));
saveas(f5,strcat(exp_ID,' processed.jpg'));

end

%% Calculate dF/F and Write to File
Delta465 = hf465./baseline465;
Delta405 = hf405./baseline405;
dFF465_M = [Ts, Delta465*100, Delta405*100]; % converts dF/F to percentages
var_array = {'Time_sec', '465_dF/F_(%)','405_dF/F_(%)',};

% Create Excel File with array data
dFF465_T = array2table(dFF465_M, 'VariableNames', var_array);
filename465 = strcat(string(exp_ID),' 465 dFF.xlsx');
writetable(dFF465_T, filename465, 'WriteVariableNames', true);

if num_sig == 4
    Delta565 = hf565./baseline565;
    dFF565_M = dFF465_M;
    dFF565_M(:,2) = Delta565*100;
    var_array{1,2} = '565_dF/F_(%)';
    
    % Create Excel File with array data
    dFF565_T = array2table(dFF565_M, 'VariableNames', var_array);
    filename565 = strcat(string(exp_ID),' 565 dFF.xlsx');
    writetable(dFF565_T, filename565, 'WriteVariableNames', true);
    
else 
    dFF565_M = [];
end

end

