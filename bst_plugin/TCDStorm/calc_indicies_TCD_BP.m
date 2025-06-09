function calc_indicies_TCD_BP(signal, RwaveLocations, signalName, result_folder_path, prefixedExtractedPattern, timePoint)

    global restStart restEnd standingStart standingEnd participant_csv_data Fs

    Fs=1000;

    % Initialize variables for segment start and end times
    segmentStart = 0;
    segmentEnd = length(signal);

    % Identify the condition ('rest' or 'standing') from the signalName
    conditionPattern = regexp(signalName, '_(rest|standing)$', 'tokens');
    if ~isempty(conditionPattern)
        condition = conditionPattern{1}{1};

        % Determine the variable names for start and end times based on the condition
        startVarName = sprintf('%sStart', condition);
        endVarName = sprintf('%sEnd', condition);

        if strcmp(condition, 'rest')
            segmentStart = restStart * Fs + 1; % Adding 1 for MATLAB indexing
            segmentEnd = restEnd * Fs;
        elseif strcmp(condition, 'standing')
            segmentStart = standingStart * Fs + 1;
            segmentEnd = standingEnd * Fs;
        end
    end

    % Adjust RwaveLocations to align with the original signal's indexing
    RwaveLocations = RwaveLocations(RwaveLocations >= segmentStart & RwaveLocations <= segmentEnd) - segmentStart + 1;

    % Check if RwaveLocations is now empty after filtering
    if isempty(RwaveLocations)
        error('Filtered RwaveLocations is empty. Check your RwaveLocations input.');
    end

    % Recalculate the number of intervals
    numIntervals = numel(RwaveLocations) - 1;

    % Initialize arrays to store calculated features
    signalAmplitude = zeros(1, numIntervals);
    signalMean = zeros(1, numIntervals);
    signalPulsatilityIndex = zeros(1, numIntervals);
    signalTimeToSystolicPeak = zeros(1, numIntervals);
    signalResistanceIndex = zeros(1, numIntervals);

    % Process the signal for all types, generic calculation logic
    for n = 1:numIntervals
        intervalStart = RwaveLocations(n);
        intervalEnd = min(RwaveLocations(n + 1), length(signal));  % Ensure intervalEnd is within bounds

        % Extract signal within the interval
        intervalSignal = signal(intervalStart:intervalEnd);

        % Perform calculations for the interval
        first40PercentEnd = min(round(0.4 * length(intervalSignal)) + intervalStart - 1, length(signal));
        [~, minIndex] = min(signal(intervalStart:first40PercentEnd));
        localMinIndex = intervalStart + minIndex - 1;

        [~, maxIndex] = max(signal(intervalStart:intervalEnd));
        localMaxIndex = intervalStart + maxIndex - 1;

        % Calculate features
        signalTimeToSystolicPeak(n) = (localMaxIndex - intervalStart) / (intervalEnd - intervalStart);
        signalAmplitude(n) = signal(localMaxIndex) - signal(localMinIndex);
        signalMean(n) = mean(intervalSignal);
        signalPulsatilityIndex(n) = (signal(localMaxIndex) - signal(localMinIndex)) / mean(intervalSignal);
        signalResistanceIndex(n) = (signal(localMaxIndex) - signal(localMinIndex)) / signal(localMaxIndex);
    end

    % Round the calculated features to 3 decimal places
    signalAmplitude = round(signalAmplitude, 3);
    signalMean = round(signalMean, 3);
    signalPulsatilityIndex = round(signalPulsatilityIndex, 3);
    signalTimeToSystolicPeak = round(signalTimeToSystolicPeak, 3);
    signalResistanceIndex = round(signalResistanceIndex, 3);

    % Dynamically adjust names based on signal type
    if contains(signalName, 'RMCA') || contains(signalName, 'LMCA')
        amplitudeUnit = '(cm/s)'; % Use cm/s for RMCA and LMCA
    else
        amplitudeUnit = '(mmHg)'; % Default to mmHg for other signals
    end
    timeUnit = '(% of cardiac cycle)';
    amplitudeName = strcat(signalName, '_amplitude', amplitudeUnit);
    meanName = strcat('mean_', signalName, amplitudeUnit);
    pulsatilityIndexName = strcat(signalName, '_pulsatilityindex');
    timeToSystolicPeakName = strcat(signalName, '_timetosystolicpeak', timeUnit);
    resistanceIndexName = strcat(signalName, '_resistanceindex');

    % Check the usability for 'rest' or 'standing' condition for LMCA and RMCA
    usabilityFlag = false;
    if strcmp(condition, 'rest') 
        if contains(signalName, 'LMCA') && isfield(participant_csv_data, 'LMCA_rest_usability') && participant_csv_data.LMCA_rest_usability == 0
            usabilityFlag = true;
        elseif contains(signalName, 'RMCA') && isfield(participant_csv_data, 'RMCA_rest_usability') && participant_csv_data.RMCA_rest_usability == 0
            usabilityFlag = true;
        end
    elseif strcmp(condition, 'standing') 
        if contains(signalName, 'LMCA') && isfield(participant_csv_data, 'LMCA_standing_usability') && participant_csv_data.LMCA_standing_usability == 0
            usabilityFlag = true;
        elseif contains(signalName, 'RMCA') && isfield(participant_csv_data, 'RMCA_standing_usability') && participant_csv_data.RMCA_standing_usability == 0
            usabilityFlag = true;
        end
    end

    % Create a table to store the median values or NaNs based on usability
    if usabilityFlag
        resultTable = table(NaN, NaN, NaN, NaN, NaN, 'VariableNames', {amplitudeName, meanName, pulsatilityIndexName, timeToSystolicPeakName, resistanceIndexName});
    else
        resultTable = table(median(signalAmplitude), median(signalMean), median(signalPulsatilityIndex), median(signalTimeToSystolicPeak), median(signalResistanceIndex), ...
            'VariableNames', {amplitudeName, meanName, pulsatilityIndexName, timeToSystolicPeakName, resistanceIndexName});
    end

    fprintf('Note: All variables reported as their median values or NaN if usability is 0.\n');

    % Save resultTable as a CSV file
    outputFileName = strcat(signalName, '_results.csv');
    writetable(resultTable, fullfile(result_folder_path, outputFileName));
    fprintf('Results saved to %s\n', fullfile(result_folder_path, outputFileName));
end
