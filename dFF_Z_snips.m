function [snips_T, start_time] = dFF_Z_snips(dFF_M, Fs, event_times, exp_ID, figtoggle, start_time)
%% Event-locked Analysis to find Maximum dF/F Values and Z-scores
% dFF_Z_snips uses a pre-defined list of event_times to generate a series
% of snips from the intact dFF trace. The function generates a heat map
% that covers all events in the session and outputs a table with the onset
% time, maximum dF/F, maximum Z-score, and peak latency for each event.
% James Maksymetz and Max Joffe April 2020

%% Inputs
%
% # dFF_M, matrix containing, times (1), Ch465 dFF (2), and Ch405 dFF (3) 
% # Fs, sampling frequency
% # Block (string) containing photometry data
% # Toggle for figures (1 = print)

%% Output
%
% snips_T, analyzed data table
%
% # Absolute times of behavioral events
% # Maximum dF/F occuring after each event
% # Maximum Z-score occuring after each event
% # Time of peak relative to behavioral event onset

% dum_vars = {'dFF', 'Fs', 'Ts', 'Block'};
% load('sub7 workspace.mat', dum_vars{:}); % sample experiment for publishing

%% Polymorphic options

if nargin == 1
        error('Not enough input variables to proceed.')
        
elseif nargin < 6
        start_time = 0;
        start_time = input('What time does session start (s)? ');    
            
        if nargin < 5
            figtoggle = 0; % figtoggle default offstart_time = 0;
            start_time = input('What time does session start (s)? ');
            
            if nargin < 4
                exp_ID = char(datetime('now')); % timestamp will replace Block if not provided
                exp_ID = exp_ID(exp_ID~=':'); % remove colons
                
                if nargin < 3
                    listFolder = uigetdir(cd, 'Select Folder with Event List');
                    cd(listFolder);
                    addpath(listFolder); 
                    event_file = uigetfile('*.csv'); % get .csv file name
                    event_times = readtable(event_file,'Delimiter',',','ReadVariableNames', false); % import .csv file
                    event_times = table2array(event_times); % convert table to array
                                 
                    if nargin == 0
                        BlockDir=uigetdir(cd, 'Select TDT Photometry Block'); % get folder name
                        cd(BlockDir); % change directory to photometry folder
                        addpath(BlockDir);
                        [Data, Fs, exp_ID] = extract_TDT_Tank(BlockDir); % extract raw data from Tank
                        [dFF_M] = extract_dFF(Data, Fs, exp_ID, figtoggle); % filter and convert data to dFF
                    end
                end
            end
        end
end

%% Set Analysis Specifications
% Customizable analysis parameters to select maximum search windows and
% baseline periods relative to behavioral events

% Set maximum length of session (seconds)
session_time = 1200;

% Minimum interval between behavioral events
min_interval = 5;

% dF/F and Z-score maximums are calculated as the maximum value within 
% peak_window (default 3 sec), after earch behavioral event onset 
peak_window = 5;

% dF/F maximum is relative to the mean during base_per_dFF (-0.5 to 0 sec)
base_per_dFF = [-0.5 0];

% Z-score baselines for standard deviation and mean are taken from
% base_per_Z (-10 to 0 sec)
base_per_Z = [-10, 0];

% Set range specifications
pre = 10; % time before behavioral event onset
post = 5; % time after behavioral event onset
range = [floor(-1*pre*Fs),ceil(post*Fs)];

% Relative location of time and dFF data in dFF_M matrtix. Change this
% value to analyze 405 or 565 channel.
time_col = 1;
dFF_col = 2;

% Check user specifications for errors
if peak_window > post 
    error('peak_window must lie within range specifications.') 
elseif base_per_Z(1) < -1*pre || base_per_Z(2) > post
    error('base_per_dFF must lie within range specifications.') 
elseif ~isa(event_times,'double')
    error('event_file must be imported as a table and converted to double.')
end

%% Check Behavioral Event List

% Ensure all behavioral events fall within session and analysis range
event_times = event_times(event_times > start_time + pre);
event_times = event_times(event_times <= start_time + session_time - post);

% Loop to ensure behavioral events are separated by min_interval
del_events = 0; % counter for deleted event values
event_times = sort(event_times);
for ii = 2:length(event_times)-1
    if event_times(ii) < event_times(ii-1) + min_interval
        event_times(ii) = 0;
        del_events = del_events + 1;
    end
end

event_times = nonzeros(event_times);

%% Split Data Matrix into Time Vector and dFF Vector
Ts = dFF_M(:,time_col); % vector containing times
dFF = dFF_M(:,dFF_col); % dFF 465 signal

%% Allocate Memory for Peri-event Snips
trials = numel(event_times);  % number of events
dFF_snips = cell(trials,1);   % cell containing a dFF vector for each event
array_ind = zeros(trials,1);  % time index corresponding to each event onset
pre_event = zeros(trials,1);  % time index leading each event onset
post_event = zeros(trials,1); % time index lagging each event onset

%% Extract Snippets of Data from dFF Trace Surrounding Each Behavioral Event
% dFF_snips corresponds to snippets of dFF values surrounding each event

for i = 1:trials
   
   % Find first time index after bout onset
   array_ind(i) = find(Ts > event_times(i),1);
 
   % Find index corresponding to pre and post stim durations
   pre_event(i) = array_ind(i) + range(1);
   post_event(i) = array_ind(i) + range(2);
   dFF_snips{i} = dFF(pre_event(i):post_event(i));
   
end

% Ensure each vector in dFF_snips is the same size
minLength = min(cellfun('prodofsize', dFF_snips));
dFF_snips = cellfun(@(x) x(1:minLength), dFF_snips, 'UniformOutput',false);

% Convert dFF_snips to a matrix
dFF_snips_mat = cell2mat(dFF_snips');
dFF_snips_mat = dFF_snips_mat';

% Make a relative time vector snippet for peri-events
peri_time = (1:length(dFF_snips_mat(1,:)))/Fs - pre;


%% Determine Maximum dFF After Event Onset

% maxd_FF_pos is maximimum dFF during each snip
% AUC_dFF is the area under the curve during peak_window (units dF/F*sec)
% max_dFF_Times is the time of each peak relative to event 

[max_dFF_pos, max_dFF_Times] = max(dFF_snips_mat(:, round(pre*Fs):round((pre+peak_window)*(Fs))),[],2); 
baseline_dFF = mean(dFF_snips_mat(:, round(((pre + base_per_dFF(1))*Fs)):round((pre + base_per_dFF(2))*Fs)),2);
max_dFF = max_dFF_pos - baseline_dFF;
max_dFF = max_dFF';
max_dFF_Times = (((max_dFF_Times + round(pre*Fs))/Fs)-pre)';
AUC_dFF = sum((dFF_snips_mat(:, round(pre*Fs):round((pre+peak_window)*Fs))-baseline_dFF),2)...
    *peak_window/1000;

%% Z-score analysis

z_snips = zeros(size(dFF_snips_mat));
for i = 1:size(dFF_snips_mat,1)
    ind = peri_time(1,:) < base_per_Z(2) & peri_time(1,:) > base_per_Z(1);
    zb = mean(dFF_snips_mat(i,ind)); % baseline period mean
    zsd = std(dFF_snips_mat(i,ind)); % baseline period stdev
    z_snips(i,:)=(dFF_snips_mat(i,:) - zb)/zsd; % Z-score per snip
end

% Determine maximum Z-score within peak_window after event onset
max_Z = max(z_snips(:, round(pre*Fs):round((pre+peak_window)*Fs)),[],2);
%base_Z = mean(z_snips(:, round((pre+base_per_dFF(1))*Fs):round((pre+base_per_dFF(2))*Fs)),2);
%net_Z = max_Z - base_Z;

%% Make a Heat Map for dF/F and Z-score Across Each Event
    
% dF/F Subplot
f2 = figure;
subplot(2,1,1)
imagesc(peri_time, 1, dFF_snips_mat); 
set(gca,'YDir','normal') % re-orders y-axis
title([exp_ID,' Struggle Bout Heat Map'],'fontsize',16)
ylabel('Trial Number','fontsize',16)
cb = colorbar;
ylabel(cb, 'dF/F (%)','fontsize',16)
axis tight;

% Z-Score Subplot
subplot(2,1,2)
imagesc(peri_time, 1, z_snips);
set(gca,'YDir','normal') % re-order y-axis
colormap(jet) % colormap otherwise defaults to parula
ylabel('Trial Number','fontsize',16)
xlabel('Seconds from struggle onset','fontsize',16)
cb = colorbar;
ylabel(cb, 'Z-Score','fontsize',16)
axis tight;

% Write Figure to file
if figtoggle == 1
savefig(strcat(exp_ID,' heat map'));
saveas(f2,strcat(exp_ID,' heat map.jpg'));
end

%% Output Maximum Values to Spreadsheet

% Create matrix with relevant values 
snips_M = [event_times, max_dFF', AUC_dFF, max_Z, max_dFF_Times'];
snips_var_names = {'Time_sec', 'Max_dFF','AUC_dFF','Max_Z_Score','Peak_Latency_sec'};

% Create table and write to Excel file
snips_T = array2table(snips_M, 'VariableNames', snips_var_names);
filename =  strcat(exp_ID,' DLC analysis.xlsx');
writetable(snips_T, filename, 'Sheet', 1);

%close all

