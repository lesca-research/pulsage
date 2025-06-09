function analyse_FingerBP(TCD_manifest_folder_path, fname, result_folder_path)

% Code is written by Hanieh Mohammadi. Version 28 March 2024.
% For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca
    

% Call the subfunction to get necessary variables from the workspace
    call_variables_from_workspace(TCD_manifest_folder_path, fname);

    % Declare global variables
    global unit_adjusted_TCD_data hasRest hasStanding restStart restEnd standingStart standingEnd RwaveLocations titles timePoint participant_csv_data Fs prefixedExtractedPattern
    
    % Check global variables for existence
    globalVars = {'unit_adjusted_TCD_data', 'hasRest', 'hasStanding', 'restStart', 'restEnd', 'standingStart', 'standingEnd', 'RwaveLocations', 'titles'};
    for varName = globalVars
        if isempty(eval(varName{1}))
            error('%s variable not found in the workspace.', varName{1});
        end
    end

    % Convert 'titles' to a cell array and find the index for 'Finger BP HCorr'
    titlesCell = cellstr(titles);
    FingerBPRowIndex = find(contains(titlesCell, 'Finger BP HCorr', 'IgnoreCase', true));
    HCUIndex = []; % Initialization for HCU index

    % Logic to handle uncorrected Finger BP and HCU data
    if isempty(FingerBPRowIndex)
        FingerBPRowIndex = find(contains(titlesCell, 'Finger BP', 'IgnoreCase', true));
        HCUIndex = find(contains(titlesCell, 'HCU', 'IgnoreCase', true));
        if isempty(FingerBPRowIndex)
            error('FingerBP data not found.');
        else
            if ~isempty(HCUIndex) % If HCU data is also found
                fprintf('Uncorrected Finger BP data found. Adding HCU data to create corrected value of Finger BP.\n');
            else
                fprintf('Using uncorrected Finger BP data. You need to height correct the FingerBP.\n');
            end
        end
    end

    % Extract the time point from the file name (e.g., "T12")
    timePointPattern = regexp(fname, 'T\d+', 'match');
    if isempty(timePointPattern)
        error('Time point not found in the file name.');
    end
    timePoint = timePointPattern{1};

    % Configuration array for periods (rest and standing)
    periods = {
        struct('flag', hasRest, 'startVar', restStart, 'endVar', restEnd, 'signalName', 'FingerBP_rest', 'filePrefix', 'FingerBP_rest'),
        struct('flag', hasStanding, 'startVar', standingStart, 'endVar', standingEnd, 'signalName', 'FingerBP_standing', 'filePrefix', 'FingerBP_standing')
    };

    % Iterating over each period configuration
    for i = 1:length(periods)
        period = periods{i};
        if period.flag % Check if period is available
            % Convert period start and end times to sample indices
            startIdx = round(period.startVar * 1000) + 1;
            endIdx = round(period.endVar * 1000);
            
            % Adjust indexing for correct data access
            signalSegment = unit_adjusted_TCD_data{FingerBPRowIndex}(startIdx:endIdx);
            if ~isempty(HCUIndex) % Add HCU data to Finger BP if HCU index is found
                signalSegment = signalSegment + unit_adjusted_TCD_data{HCUIndex}(startIdx:endIdx);
            end

            % Process the signal segment
            calc_indicies_TCD_BP(signalSegment, RwaveLocations, period.signalName, result_folder_path, prefixedExtractedPattern, timePoint);
            % File renaming logic
            originalFileName = sprintf('%s_results.csv', period.signalName);
            newFileName = sprintf('%s_%s_%s_results.csv', prefixedExtractedPattern, period.filePrefix, timePoint);
            originalFilePath = fullfile(result_folder_path, originalFileName);
            newFilePath = fullfile(result_folder_path, newFileName);
    
            % Rename the file if it exists
            if exist(originalFilePath, 'file')
                movefile(originalFilePath, newFilePath);
                fprintf('File renamed to: %s\n', newFileName);
            else
                fprintf('Original file for %s does not exist, cannot rename.\n', period.signalName);
            end
        end
    end
end
