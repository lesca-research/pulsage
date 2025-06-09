function ECG_RwaveLocations(TCD_manifest_folder_path, fname)

% Code is written by Hanieh Mohammadi. Version 28 March 2024.
% For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca

% Declare global variables to share ECG data and R-wave locations across functions.
global unit_adjusted_TCD_data RwaveLocations participant_csv_data;

%% Code is written by Hanieh Mohammadi Version 19 Feb 2024
%% For any question please contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca

% Retrieve ECG analysis parameters from a specified CSV file.
call_variables_from_workspace(TCD_manifest_folder_path, fname);

% Extract relevant parameters from the CSV data for ECG analysis.
ECG_leadNum = participant_csv_data.selected_ECG_lead_num; % The ECG lead number to analyze.
ECG_MinPeakHeight = participant_csv_data.ECG_MinPeakHeight; % Minimum height of R-wave peaks.
ECG_minPeakDistance = participant_csv_data.ECG_minPeakDistance; % Minimum distance between R-wave peaks.
negateECG = participant_csv_data.negate_selected_ecg; % Flag to check if ECG data should be negated.

% Check if ECG_leadNum is zero or NA (indicating not applicable)
if isequal(ECG_leadNum, 0) || isequal(lower(ECG_leadNum), 'na')
    fprintf('ECG is not usable. Please tag the BP and LMCA / or LMCA manually.\n');
    return; % Exit the function since ECG analysis cannot proceed.
end

% Ensure the ECG data variable is not empty before proceeding.
if isempty(unit_adjusted_TCD_data)
    fprintf('The ECG file (unit_adjusted_tcd_data) is not in the workspace or is empty.\n');
    return; 
end

% Validate the ECG lead number and extract the corresponding ECG data.
if ECG_leadNum < 1 || ECG_leadNum > size(unit_adjusted_TCD_data, 1)
    fprintf(['Invalid ECG lead number: ', num2str(ECG_leadNum)]); % Inform the user if the lead number is invalid.
    return; % Exit the function if the lead number is invalid.
end

original_ecg = unit_adjusted_TCD_data{ECG_leadNum, 1}; % Extract ECG data for the specified lead.
disp(['ECG data for lead number ', num2str(ECG_leadNum)]); % Display the lead number being analyzed.

% Define the ECG sampling frequency and create a time vector for plotting.
Fs = 1000; % Sampling frequency in Hz.
t = (0:length(original_ecg)-1) / Fs; % Time vector is defined here, ensuring it's available for plotting.

% Negate the ECG data for peak detection if required.
ecg_for_detection = original_ecg;
if negateECG == 1
    ecg_for_detection = -original_ecg;
end

% Use the `findpeaks` function to identify R-wave peaks on the potentially negated data.
[~, RwaveLocations] = findpeaks(ecg_for_detection, 'MinPeakHeight', ECG_MinPeakHeight, 'MinPeakDistance', ECG_minPeakDistance*Fs);

% Plot the original ECG signal and mark R-wave locations.
figure; 
plot(t, original_ecg); % Use original ECG data for plotting, ensured 't' is defined before this line
hold on; 
plot(t(RwaveLocations), original_ecg(RwaveLocations), 'ro', 'MarkerFaceColor', 'r');
title('ECG Signal with R-wave Peaks'); 
xlabel('Time (s)'); 
ylabel('ECG'); 
legend('ECG Signal', 'R-wave Peaks');

end
