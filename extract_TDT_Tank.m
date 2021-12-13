function [raw_data, Fs, exp_ID] = extract_TDT_Tank(BlockDir)
%% Extract Data from TDT Tank 
% Process_TDT_Tank imports data from a TDT Tank folder using the TDTbin2mat
% function. Data is then down-sampled (ten-fold, by default) and exported
% as an array with first row corresponding to timestamps and each following
% vector corresponding to a fluorescence channel. 
% James Maksymetz and Max Joffe April 2020

%% Input
% Process_TDT_Tank can accept BlockDir, the directory of the TDT Tank as an
% input. If no input is applied, Process_TDT_Tank calls for the user to
% select a directory with a user interface

% dummy Tank directory'C:\191025_SST_7-191122-123444\'

%% Outputs
%
% # raw_data, array containing continuous data (row 1:time, 2: 465, 3:405, 4:565)
% # Fs, final sampling frequency

%% Load Data Using TDTbin2mat
if nargin == 0
BlockDir=uigetdir(cd, 'Select TDT Photometry Block'); % get folder name
end
cd(BlockDir); % change directory to photometry folder

data = TDTbin2mat(BlockDir); 

exp_ID = data.info.blockname;

%% Extract Information from Data Structure

channel_names = fieldnames(data.streams); % Assign channel names to a cell array

for ii = 1:length(fieldnames(data.streams))
   if strfind(channel_names{ii},'405') > 0 % Strings to substitute in data.streams.xxxx.data
       BLUE = channel_names{ii};
   elseif strfind(channel_names{ii},'465') > 0
       GREEN = channel_names{ii};
   elseif strfind(channel_names{ii},'565') > 0 
       RED = channel_names{ii};
   end    
end

% Downsample data to speed up processing
N = 10; % multiplicative for downsampling
data.streams.(GREEN).data = arrayfun(@(i)...
    mean(data.streams.(GREEN).data(i:i+N-1)),...
    1:N:length(data.streams.(GREEN).data)-N+1);
data.streams.(BLUE).data = arrayfun(@(i)...
    mean(data.streams.(BLUE).data(i:i+N-1)),...
    1:N:length(data.streams.(BLUE).data)-N+1);

% Create variables for each channel
Ch465=data.streams.(GREEN).data; % GCaMP
Ch405=data.streams.(BLUE).data; % isosbestic control

% Determine sampling frequency and timestamps
Ts = ((1:numel(data.streams.(GREEN).data(1,:))) / (data.streams.(GREEN).fs/N))'; % get Ts for samples based on Fs
StartTime = 100; % set the starting sample(recommend eliminating a few seconds for photoreceiver/LED rise time)
EndTime = length(Ch465) - 100; % set the ending sample (again, eliminate some)
Fs=data.streams.(GREEN).fs/N; % sampling frequency
Ts=Ts(StartTime:EndTime); % eliminate timestamps before starting sample and after ending.

% Trim data files
Ch465 = Ch465(StartTime:EndTime);
Ch405 = Ch405(StartTime:EndTime);

% Concatenate data into one array to export
raw_data = [Ts'; Ch465; Ch405];

%% Extract Red Channel if It Exists
if numel(channel_names) == 4
    data.streams.(RED).data = arrayfun(@(i)...
        mean(data.streams.(RED).data(i:i+N-1)),...
        1:N:length(data.streams.(RED).data)-N+1);

    Ch565 = data.streams.(RED).data; % RCaMP or tomato control
    Ch565 = Ch565(StartTime:EndTime);
    raw_data = [raw_data; Ch565];
end

close all

end
