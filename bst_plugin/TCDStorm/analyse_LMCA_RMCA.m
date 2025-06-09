function analyse_LMCA_RMCA(TCD_manifest_folder_path, fname, result_folder_path)
    
    % Code is written by Hanieh Mohammadi. Version 28 March 2024.
    % For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca

    % Declare global variables
    global unit_adjusted_TCD_data hasRest hasStanding restStart restEnd standingStart standingEnd RwaveLocations titles prefixedExtractedPattern Fs

    % Call the subfunction to get necessary variables from the workspace
    call_variables_from_workspace(TCD_manifest_folder_path, fname);

    % Check for the existence of global variables
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
        error(['The following required variable(s) are not found in the workspace: ' errorStr '.']);
    end

    disp('All required global variables are found in the workspace.');


    titlesCell = cellstr(titles);
    LMCARowIndex = find(contains(titlesCell, 'TCD L (red)', 'IgnoreCase', true)); % Identifying LMCA data
    RMCARowIndex = find(contains(titlesCell, 'TCD R (green)', 'IgnoreCase', true)); % Identifying RMCA data
    
    % Check if specific data types are available
    if isempty(LMCARowIndex) && isempty(RMCARowIndex)
        error('Neither LMCA nor RMCA data found.');
    end

    % Extract the time point from fname (e.g., "T12")
    timePointPattern = regexp(fname, 'T\d+', 'match');
    if isempty(timePointPattern)
        error('Time point not found in the file name.');
    end
    timePoint = timePointPattern{1};

    % Dynamically create period configurations based on available data
    periods = {};
    if hasRest
        if ~isempty(LMCARowIndex)
            periods{end+1} = struct('startVar', 'restStart', 'endVar', 'restEnd', 'signalName', 'LMCA_rest', 'filePrefix', 'LMCA_rest', 'rowIndex', LMCARowIndex);
        end
        if ~isempty(RMCARowIndex)
            periods{end+1} = struct('startVar', 'restStart', 'endVar', 'restEnd', 'signalName', 'RMCA_rest', 'filePrefix', 'RMCA_rest', 'rowIndex', RMCARowIndex);
        end
    end
    if hasStanding
        if ~isempty(LMCARowIndex)
            periods{end+1} = struct('startVar', 'standingStart', 'endVar', 'standingEnd', 'signalName', 'LMCA_standing', 'filePrefix', 'LMCA_standing', 'rowIndex', LMCARowIndex);
        end
        if ~isempty(RMCARowIndex)
            periods{end+1} = struct('startVar', 'standingStart', 'endVar', 'standingEnd', 'signalName', 'RMCA_standing', 'filePrefix', 'RMCA_standing', 'rowIndex', RMCARowIndex);
        end
    end

    % Process each available period for LMCA and RMCA
    for i = 1:length(periods)
        period = periods{i};
        % Convert period to sample indices
        % Directly use global variables based on the condition
            if strcmp(period.startVar, 'restStart')
                startIdx = round(restStart * 1000) + 1;
            elseif strcmp(period.startVar, 'standingStart')
                startIdx = round(standingStart * 1000) + 1;
            end
            
            if strcmp(period.endVar, 'restEnd')
                endIdx = round(restEnd * 1000);
            elseif strcmp(period.endVar, 'standingEnd')
                endIdx = round(standingEnd * 1000);
            end
        signalSegment = unit_adjusted_TCD_data{period.rowIndex}(startIdx:endIdx);

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
