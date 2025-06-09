close all;
clear all;
clc;

%% Code information
% Code is written by Hanieh Mohammadi
% Version 21 March 2024
% For any questions, please contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca

%% Setup file paths and names
fname = 'RC088_T0_rest_stand_TCD_RMCA_LMCA.mat';
path_to_LabChart_data = '/Users/haniehm/Documents/ACTIONcrNov1st2023/ACTION_TCD/TCDStorm/TroubleShooting_Recardio/';
csv_folder_path = '/Users/haniehm/Documents/ACTIONcrNov1st2023/ACTION_TCD/TCDStorm/TroubleShooting_Recardio/';
csv_filename = fullfile(csv_folder_path, 'TCD_manifest_main_test_modified.csv');

%% Read CSV for ECG parameters using readtable directly
opts = detectImportOptions(csv_filename, 'NumHeaderLines', 0);
T = readtable(csv_filename, opts);

% Remove the .mat extension from fname for matching purposes
fnameWithoutExt = fname(1:end-4);
subjectRowIndex = find(strcmp(T.datafile_orgName, fnameWithoutExt));

%% Extract necessary variables
selected_ECG_lead_num = T.selected_ECG_lead_num(subjectRowIndex); % Adjusted to use the subject's row
ECG_MinPeakHeight = T.ECG_MinPeakHeight(subjectRowIndex);
ECG_minPeakDistance = T.ECG_minPeakDistance(subjectRowIndex);
negate_selected_ecg = T.negate_selected_ecg(subjectRowIndex); % Assuming this exists in your CSV

% Check for missing or NaN values
if isempty(selected_ECG_lead_num) || isnan(selected_ECG_lead_num) || ...
   isempty(ECG_MinPeakHeight) || isnan(ECG_MinPeakHeight) || ...
   isempty(ECG_minPeakDistance) || isnan(ECG_minPeakDistance) || ...
   isempty(negate_selected_ecg) || isnan(negate_selected_ecg)
    error('One or more required variables are missing or NaN. The plot could not be completed.');
end

%% Load LabChart data
load(fullfile(path_to_LabChart_data, fname));

if ~exist('data','var')
    error('No data found. Ensure the .mat file contains "data" variable and is exported correctly.');
end

%% Extract and preprocess the specified ECG lead data
ECG_data = data(datastart(selected_ECG_lead_num,1):dataend(selected_ECG_lead_num,1));
ptime = (0:length(ECG_data)-1) / samplerate(selected_ECG_lead_num,1);

% Temporarily negate ECG data for peak detection if required
ECG_data_for_detection = ECG_data;
if negate_selected_ecg == 1
    ECG_data_for_detection = -ECG_data;
end

%% Detect R-waves
Fs = samplerate(selected_ECG_lead_num,1); % Sampling frequency
[~, RwaveLocations] = findpeaks(ECG_data_for_detection, 'MinPeakHeight', ECG_MinPeakHeight, 'MinPeakDistance', ECG_minPeakDistance*Fs);

%% Plot the ECG data and R-wave locations
figure;
plot(ptime, ECG_data);
hold on;
plot(RwaveLocations/Fs, ECG_data(RwaveLocations), 'r.', 'MarkerSize', 10); % R-wave locations as red dots
title(sprintf('ECG Lead %d with R-wave locations', selected_ECG_lead_num));
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 max(ptime)]);

%% Annotate and Print context on ECG Plot
% Assuming 'com' and 'comtext' are variables within your .mat file related to event timestamps and annotations
if exist('com','var') && exist('comtext','var')
    fprintf('\nEvent timestamps and annotations:\n');
    for i = 1:size(com, 1)
        comIndex = com(i,5); % Index for comtext
        eventText = strtrim(comtext(comIndex,:)); % Actual text
        eventTextEscaped = strrep(eventText, '_', '\_'); % Escape underscores
        eventTimeInSeconds = com(i,3) / Fs;
        fprintf('%s occurred at %.3f seconds\n', eventTextEscaped, eventTimeInSeconds);
        
        % Annotate on plot if within bounds
        if eventTimeInSeconds <= max(ptime)
            dataIndex = round(eventTimeInSeconds * Fs);
            if dataIndex <= length(ECG_data)
                text(eventTimeInSeconds,ECG_data(dataIndex), eventTextEscaped, 'Rotation', 90, 'BackgroundColor', 'white', 'VerticalAlignment', 'bottom');
            end
        end
    end
else
    fprintf('Event information is missing or incomplete in the provided data.\n');
end

