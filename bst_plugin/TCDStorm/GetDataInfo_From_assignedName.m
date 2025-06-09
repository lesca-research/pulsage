function GetDataInfo_From_assignedName(TCD_manifest_folder_path, fname)
  
% Code is written by Hanieh Mohammadi. Version 28 March 2024.
% For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca
 

% Declare global variables
    global hasRest restStart restEnd hasStanding standingStart standingEnd Fs
    
     Fs=1000; % hz
    % Define the name of the file
    csv_fileName = 'RC_TCD_manifest.csv';
    
    % Use fullfile to construct the full file path
    fullFilePath = fullfile(TCD_manifest_folder_path, csv_fileName);
    
    % Read the CSV file into a table
    csv_data_Table = readtable(fullFilePath);
    
    % Find the row with the input filename
    participantRowIndex = find(strcmp(csv_data_Table.datafile_orgName, fname));

    if isempty(participantRowIndex)
        error('The specified filename was not found in the CSV files list.');
    end

    % Extract the assigned name from the found row
    assignedName = csv_data_Table.datafile_assignedName{participantRowIndex};
    
    % Check for "rest" and "standing" in the assigned name
    hasRest = contains(assignedName, 'rest');
    hasStanding = contains(assignedName, 'standing');
    
    % Define usability check function for conciseness
    function usable = isUsable(index, types)
        usable = any(csv_data_Table{index, types} == 1);
    end
    
    % Initialize output variables
    restStart = 0;
    restEnd = 0;
    standingStart = 0;
    standingEnd = 0;
    
    % Extract corresponding start and end times based on the assigned name
    if hasRest && isUsable(participantRowIndex, {'ECG_lead1_rest_usability', 'ECG_lead2_rest_usability', 'ECG_lead3_rest_usability'})
        restStart = csv_data_Table.rest_start(participantRowIndex);
        restEnd = csv_data_Table.rest_end(participantRowIndex);
    else
        fprintf('Rest period data is not available or usable for "%s".\n', assignedName);
    end

    if hasStanding && isUsable(participantRowIndex, {'ECG_lead1_standing_usability', 'ECG_lead2_standing_usability', 'ECG_lead3_standing_usability'})
        standingStart = csv_data_Table.standing_start(participantRowIndex);
        standingEnd = csv_data_Table.standing_end(participantRowIndex);
    else
        fprintf('Standing period data is not available or usable for "%s".\n', assignedName);
    end
    
end
