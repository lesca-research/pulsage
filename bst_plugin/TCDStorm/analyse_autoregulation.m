
function analyse_autoregulation(TCD_manifest_folder_path, fname, result_folder_path, windowSizeInSeconds)

% Code is written by Hanieh Mohammadi. Version 28 March 2024.
% For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca

%% moving_average_window_in_minute could be 0 minute or 10 second or 60 sec minute or two minute

% Call the subfunction to get necessary variables from the workspace
call_variables_from_workspace(TCD_manifest_folder_path, fname);

% Declare global variables
global unit_adjusted_TCD_data hasRest hasStanding restStart restEnd standingStart standingEnd RwaveLocations titles timePoint participant_csv_data prefixedExtractedPattern windowSizeInSeconds

% Check global variables for existence
globalVars = {'unit_adjusted_TCD_data', 'hasRest', 'hasStanding', 'restStart', 'restEnd', 'standingStart', 'standingEnd', 'RwaveLocations', 'titles'};
missingVars = {};
for varName = globalVars
    if isempty(eval(varName{1}))
        missingVars{end+1} = varName{1}; % Collect names of missing variables
    end
end

% Display error for missing variables
if ~isempty(missingVars)
    errorStr = strjoin(missingVars, ', ');
    error(['The following equired variables for Autoregulation calculation are not found in the workspace: ' errorStr '.']);
else
    disp('All required global variables required for Autoregulation calculation are found in the workspace.');
end

% Convert 'titles' to a cell array
titlesCell = cellstr(titles);
FingerBPRowIndex = find(contains(titlesCell, 'Finger BP', 'IgnoreCase', true));
HCUIndex = find(contains(titlesCell, 'HCU', 'IgnoreCase', true)); % For FingerBP correction
LMCARowIndex = find(contains(titlesCell, 'TCD L (red)', 'IgnoreCase', true));
RMCARowIndex = find(contains(titlesCell, 'TCD R (green)', 'IgnoreCase', true));

% Extract the time point from the file name
timePointPattern = regexp(fname, 'T\d+', 'match');
if isempty(timePointPattern)
    error('Time point not found in the file name.');
end
timePoint = timePointPattern{1};

% Check data availability and report (function is below this code)
checkDataAvailability(FingerBPRowIndex, HCUIndex, LMCARowIndex, RMCARowIndex);

% Configuration array for periods (rest and standing) for all data types
periods = configurePeriods(hasRest, hasStanding, restStart, restEnd, standingStart, standingEnd, FingerBPRowIndex, HCUIndex, LMCARowIndex, RMCARowIndex);

% Process each available period for FingerBP, LMCA, and RMCA
processPeriods(periods, unit_adjusted_TCD_data, RwaveLocations, result_folder_path, prefixedExtractedPattern);
end

function checkDataAvailability(FingerBPRowIndex, HCUIndex, LMCARowIndex, RMCARowIndex)
if ~isempty(FingerBPRowIndex)
    fprintf('FingerBP data found.\n');
    if ~isempty(HCUIndex)
        fprintf('HCU data found for FingerBP correction.\n');
    end
else
    fprintf('FingerBP data not found.\n');
end

if ~isempty(LMCARowIndex)
    fprintf('LMCA data found.\n');
else
    fprintf('LMCA data not found.\n');
end

if ~isempty(RMCARowIndex)
    fprintf('RMCA data found.\n');
else
    fprintf('RMCA data not found.\n');
end
end

function periods = configurePeriods(hasRest, hasStanding, restStart, restEnd, standingStart, standingEnd, FingerBPRowIndex, HCUIndex, LMCARowIndex, RMCARowIndex)
% Initialize periods as an empty struct array
periods = struct('periodType', {}, 'startIndex', {}, 'endIndex', {}, 'dataIndices', {});

% Counter for adding new period configurations
periodCounter = 0;

if hasRest
    % Increment counter
    periodCounter = periodCounter + 1;

    % Ensure indices are single values or NaN, and configure rest period
    periods(periodCounter).periodType = 'Rest';
    periods(periodCounter).startIndex = restStart;
    periods(periodCounter).endIndex = restEnd;
    periods(periodCounter).dataIndices = [...
        ensureSingleValue(FingerBPRowIndex), ...
        ensureSingleValue(HCUIndex), ...
        ensureSingleValue(LMCARowIndex), ...
        ensureSingleValue(RMCARowIndex)];
end

if hasStanding
    % Increment counter
    periodCounter = periodCounter + 1;

    % Ensure indices are single values or NaN, and configure standing period
    periods(periodCounter).periodType = 'Standing';
    periods(periodCounter).startIndex = standingStart;
    periods(periodCounter).endIndex = standingEnd;
    periods(periodCounter).dataIndices = [...
        ensureSingleValue(FingerBPRowIndex), ...
        ensureSingleValue(HCUIndex), ...
        ensureSingleValue(LMCARowIndex), ...
        ensureSingleValue(RMCARowIndex)];
end
end

function value = ensureSingleValue(index)
% Ensure the index is either a single value or NaN
if isempty(index)
    value = NaN; % Use NaN for missing indices
else
    value = index(1); % Use the first value if there are multiple (or if it's already a single value)
end
end

function processPeriods(periods, unit_adjusted_TCD_data, RwaveLocations, result_folder_path, prefixedExtractedPattern, windowSizeInSeconds)
    
global Fs participant_csv_data windowSizeInSeconds

    % Initialize variables to hold autoregulation indices for rest and standing as NaN
    AutoReg_static_nMxa_LMCA = []; % For Rest periods
    AutoReg_static_nMxa_RMCA = []; % For Rest periods
    AutoReg_static_nMxa_LMCA_Standing = []; % For Standing periods
    AutoReg_static_nMxa_RMCA_Standing = []; % For Standing periods

    % Design a bandpass filter for the desired frequency range (0.005 – 0.05 Hz)
    [b, a] = butter(2, [0.005 0.05] / (Fs/2), 'bandpass');

    % Loop through each period to process
    for i = 1:length(periods)
        periodType = periods(i).periodType;
        if strcmp(periodType, 'Rest') || strcmp(periodType, 'Standing')
            fingerBPIndex = periods(i).dataIndices(1);
            LMCAIndex = periods(i).dataIndices(3);
            RMCAIndex = periods(i).dataIndices(4);

            % Adjust start and end indices by sampling frequency and ensure startIndex is at least 1
            startIndex = max(1, periods(i).startIndex * Fs);
            endIndex = periods(i).endIndex * Fs;

            % Extract the relevant data range from these vectors/matrices for the current period
            fingerBPData = unit_adjusted_TCD_data{fingerBPIndex}(startIndex:endIndex);
            LMCAData = unit_adjusted_TCD_data{LMCAIndex}(startIndex:endIndex);
            RMCAData = unit_adjusted_TCD_data{RMCAIndex}(startIndex:endIndex);

            % Calculate correlations, considering window size
            correlationsLMCA = [];
            correlationsRMCA = [];
            if windowSizeInSeconds > 0
                windowSize = windowSizeInSeconds * Fs;
                for wStart = 1:windowSize:length(fingerBPData) - windowSize + 1
                    wEnd = wStart + windowSize - 1;
                    correlationsLMCA(end+1) = corr(fingerBPData(wStart:wEnd)', LMCAData(wStart:wEnd)');
                    correlationsRMCA(end+1) = corr(fingerBPData(wStart:wEnd)', RMCAData(wStart:wEnd)');
                end
            else
                correlationsLMCA = corr(fingerBPData', LMCAData');
                correlationsRMCA = corr(fingerBPData', RMCAData');
            end

            % Store or process the correlations based on the period type
            if strcmp(periodType, 'Rest')
                AutoReg_static_nMxa_LMCA(end+1) = median(correlationsLMCA);
                AutoReg_static_nMxa_RMCA(end+1) = median(correlationsRMCA);
            elseif strcmp(periodType, 'Standing')
                AutoReg_static_nMxa_LMCA_Standing(end+1) = median(correlationsLMCA);
                AutoReg_static_nMxa_RMCA_Standing(end+1) = median(correlationsRMCA);
            end
        end
    end

    % Calculate the dynamic indices (Standing - Rest) / Rest
    AutoReg_dynamic_nMxa_LMCA = (AutoReg_static_nMxa_LMCA_Standing - AutoReg_static_nMxa_LMCA) ./ AutoReg_static_nMxa_LMCA;
    AutoReg_dynamic_nMxa_RMCA = (AutoReg_static_nMxa_RMCA_Standing - AutoReg_static_nMxa_RMCA) ./ AutoReg_static_nMxa_RMCA;

    % Prepare to save results, including dynamic indices
    subjectCode = participant_csv_data.SubjectCode;
    timePoint = participant_csv_data.timePoint;
    filename = fullfile(result_folder_path, sprintf('%s_AutoReg_%s_results.csv', subjectCode, timePoint));

    % Create the table with all variables
    resultsTable = table(...
        num2cell(AutoReg_static_nMxa_LMCA), ...
        num2cell(AutoReg_static_nMxa_RMCA), ...
        num2cell(AutoReg_dynamic_nMxa_LMCA), ...
        num2cell(AutoReg_dynamic_nMxa_RMCA), ...
        'VariableNames', {...
        'AutoReg_static_nMxa_LMCA_Rest', ...
        'AutoReg_static_nMxa_RMCA_Rest', ...
        'AutoReg_dynamic_nMxa_LMCA', ...
        'AutoReg_dynamic_nMxa_RMCA'});

    % Write the table to a CSV file
    writetable(resultsTable, filename);

    % Report that the file was saved
    fprintf('Results with Static and Dynamic Autoreg indices based on rest and standing data saved to %s\n', filename);
end


% to do, we can add rate of regulation as an index, which is rate of change
% in CBV normalized by arterial blood pressure measured



