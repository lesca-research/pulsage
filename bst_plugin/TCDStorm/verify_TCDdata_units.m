function  verify_TCDdata_units(TCD_manifest_folder_path, fname)
   
% Code is written by Hanieh Mohammadi. Version 28 March 2024.
% For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca


    call_variables_from_workspace(TCD_manifest_folder_path, fname);

    global titles unittextmap unittext pdatas unit_adjusted_TCD_data

    % variables of interest
    searchTitles = {'Finger BP', 'TCD L', 'TCD R', 'HCU'};

    % Convert the character array 'titles' to a cell array of strings
    titlesCell = cellstr(titles);

    % Initialize an array to store row indices
    allRowIndices = [];

    % Iterate through each search title
    for k = 1:numel(searchTitles)
        % Remove all spaces from both titles and the current searchTitle
        titlesNoSpaces = cellfun(@(x) strrep(x, ' ', ''), titlesCell, 'UniformOutput', false);
        searchTitleNoSpaces = strrep(searchTitles{k}, ' ', '');

        % Find the row indices for the specified title
        rowIndices = [];
        for i = 1:length(titlesNoSpaces)
            if contains(titlesNoSpaces{i}, searchTitleNoSpaces, 'IgnoreCase', true)
                rowIndices = [rowIndices, i];
            end
        end

        % Save the row indices
        allRowIndices = [allRowIndices, rowIndices];

        % Check and convert units based on the search title
        for i = 1:length(rowIndices)
            rowIndex = rowIndices(i);

            % Fetch the unit code from unittextmap for the specified row
            unitCodeFromMap = unittextmap(rowIndex);

            % Map unit code to the corresponding unit in unittext
            unitsFromMap = unittext(unitCodeFromMap, :);

            % Check and perform conversion based on the search title
            if strcmpi(searchTitleNoSpaces, 'FingerBP') && ~strcmpi(unitsFromMap, 'mmHg')
                fprintf('Converting Volt to mmHg for FingerBP in row %d.\n', rowIndex);
                pdatas{rowIndex, 1} = 100 * pdatas{rowIndex, 1}; % Update the original data
                unitsFromMap = 'mmHg'; % Update the units

                % Add the new condition for HCU here
             elseif strcmpi(searchTitleNoSpaces, 'HCU') && ~strcmpi(unitsFromMap, 'mmHg')
                fprintf('Converting Volt to mmHg for HCU in row %d.\n', rowIndex);
                pdatas{rowIndex, 1} = 100 * pdatas{rowIndex, 1}; % Assume the same conversion factor for HCU as for FingerBP
                unitsFromMap = 'mmHg'; % Update the units

                % TO DO: add adjustment of finapres

            elseif (strcmpi(searchTitleNoSpaces, 'TCDL') || strcmpi(searchTitleNoSpaces, 'TCDR')) && ~strcmpi(unitsFromMap, 'cm/s')
                fprintf('Converting Volt to cm/s units for %s in row %d.\n', searchTitles{k}, rowIndex);
                pdatas{rowIndex, 1} = (pdatas{rowIndex, 1} - 0) / 0.97 * 100; % Update the original data
                unitsFromMap = 'cm/s'; % Update the units
            end
        end
    end

    % Assuming you want to return the updated pdatas
    unit_adjusted_TCD_data = pdatas;
end
