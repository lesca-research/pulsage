function  call_variables_from_workspace(csv_folder_path, fname)
  
global participant_csv_filename participant_csv_data assignedName prefixedExtractedPattern % Declare as global
    
    TCD_manifest_fileName = 'RC_TCD_manifest.csv';

    TCD_csv_FilePath = fullfile(csv_folder_path, TCD_manifest_fileName);
    csv_dataTable = readtable(TCD_csv_FilePath);
    participantRowIndex = find(strcmp(csv_dataTable.datafile_orgName, fname));
    
    if isempty(participantRowIndex)
        error('Participant filename not found in the CSV data.');
    end
    
    assignedName = csv_dataTable.datafile_assignedName{participantRowIndex};
    assignedName_pattern ='RC(\d+)_TCD_T(\d+)';


    tokens = regexp(assignedName, assignedName_pattern, 'tokens');
    
    if isempty(tokens)
        error('The assigned name does not contain the expected pattern.');
    end
    
    extractedNumber = tokens{1}{1};
    participant_csv_filename = ['csv_data_ACR_' extractedNumber];
    
    if evalin('base', ['exist(''' participant_csv_filename ''', ''var'')'])
        participant_csv_data = evalin('base', participant_csv_filename);
    else
        participant_csv_data = table2struct(csv_dataTable(participantRowIndex, :));
   
    end

    % Use the returned 'assignedName' to construct the prefixedExtractedPattern and other variables if necessary
    pattern = 'RC(\d+)_TCD_T(\d+)'; % Updated pattern to match the original function's needs
    tokens = regexp(assignedName, pattern, 'tokens');
    if isempty(tokens)
        error('The pattern RCx-Tx_... could not be found in assignedName.');
    end
    % Correctly extracting the first and second part of the pattern
    rcNumber = tokens{1}{1}; % Extracting the RC number
    tNumber = tokens{1}{2}; % Extracting the T number directly
    
    % Constructing the string with actual numbers
    prefixedExtractedPattern = strcat("RC", rcNumber, "_TCD_T", tNumber);
    
    % If you need it as a string object rather than char array (depending on your MATLAB version)
    prefixedExtractedPattern = string(prefixedExtractedPattern);
    
    % No need for assignin since we're using global variables
end
