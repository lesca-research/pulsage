function fetch_TCD_data(path_to_LabChart_data, TCD_manifest_folder_path, fname)
    
% Code is written by Hanieh Mohammadi. Version 28 March 2024.
% For any question, contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca


% Declare global variables
    global pdatas titles unittextmap unittext Fs oxySoftTrigger

    Fs = 1000;

    % Load the specified .mat file
    load(fullfile(path_to_LabChart_data, fname));

    % Check if 'data' variable exists in the loaded file
    if ~exist('data', 'var')
        error('No data! Select a mat file that contains data and was created with Export Matlab 3.0 or later (LabChart for Windows 7.2 or later)');
        return;
    end

    [numchannels, numblocks] = size(datastart);
    numblocks = 1; % Modification if only using one block

    % Set up figure, size, position
    figure(1);
    clf;
    FigName = ['LabChart - All Channels, all Blocks of ' fname];
    scrsz = get(0, 'ScreenSize');
    set(gcf, 'Name', FigName, 'Position', scrsz([3 4 3 4]) .* [1 1 6 6] / 8);

    % Initialize data arrays
    pdatas = cell(numchannels, numblocks);
    ptimes = cell(numchannels, numblocks);

    % Processing each channel and block
    for ch = 1:numchannels
        for bl = 1:numblocks
            if (datastart(ch, bl) ~= -1) % Check for non-empty blocks
                % Data extraction and scaling
                pdatas{ch, bl} = data(datastart(ch, bl):dataend(ch, bl));
                if exist('scaleunits', 'var')
                    pdatas{ch, bl} = (pdatas{ch, bl} + scaleoffset(ch, bl)) .* scaleunits(ch, bl);
                end
                ptimes{ch, bl} = [0:size(pdatas{ch, bl}, 2) - 1] / samplerate(ch, bl);

                % Plotting data
                subplot(numchannels, numblocks, (ch - 1) * numblocks + bl);
                plot(ptimes{ch, bl}, pdatas{ch, bl}), hold on;

                % Further processing for each block
                pdata = pdatas{ch, bl};
                ptime = ptimes{ch, bl};

                % Setting titles and labels
                if bl == 1
                    title(titles(ch, :));
                end
                xlabel('Time (s)');
                if length(ptime) ~= 1
                    xlim([0 max(ptime)]);
                end

                % Setting y-axis labels / units
                if (unittextmap(ch, bl) ~= -1)
                    unit = unittext(unittextmap(ch, bl), :);
                    ylabel(unit);
                end

                % Setting plot limits
                if exist('scaleunits', 'var')
                    pmax = (rangemax(ch, bl) + scaleoffset(ch, bl)) .* scaleunits(ch, bl);
                    pmin = (rangemin(ch, bl) + scaleoffset(ch, bl)) .* scaleunits(ch, bl);
                else
                    pmax = max(pdata);
                    pmin = min(pdata);
                end
                ylim([pmin pmax]);

                % Plotting comments if they exist
                if exist('com', 'var')
                    temp = find(com(:, 2) == bl & (com(:, 1) == ch | com(:, 1) == -1));
                    comtickpos = com(:, 3);
                    comtextmap = com(:, 5);
                    for m = 1:size(temp)
                        temp2 = temp(m);
                        x = round(comtickpos(temp2) * (samplerate(ch, bl) / tickrate(bl)));
                        t = x / samplerate(ch, bl);
                        plot([t t], [pmin pmax], 'r:');
                        
                        % Modify the comtext string to replace underscores with '\_' for LaTeX interpretation
                        modifiedText = strrep(comtext(comtextmap(temp2), :), '_', '\_');
                        
                        % Use the text function with the LaTeX interpreter to display the modified text
                        text(t, (pmin + 0.1 * (pmax - pmin)), modifiedText, 'Rotation', 90, 'Interpreter', 'latex');
                        
                        % Check if the modifiedText is 'oxysoft sync' and store its time globally
                        if strcmp(strtrim(modifiedText), 'oxysoft sync') % Corrected for LaTeX interpretation
                           oxySoftTrigger = t; % Storing the time associated with 'oxysoft sync'
                            % Here, t is the time in seconds
                        end
                    end
                end

            end
        end
    end
    
    % After processing all channels and blocks, you can optionally print or
    % take further actions with the global variable 'oxySoftTrigger' if it was set.
    if exist('oxySoftTrigger', 'var')
        disp(['OxySoft Trigger Time: ', num2str(oxySoftTrigger), ' seconds']);
    else
        disp('No OxySoft Trigger Time Found.');
    end

    % Initialize a variable to store the maximum time
    maxTime = 0;
    
    % Iterate through each channel and block to find the maximum time value
    for ch = 1:numchannels
        for bl = 1:numblocks
            % Check if the current cell in ptimes is not empty
            if ~isempty(ptimes{ch, bl})
                % Update maxTime if the last element of the current time array is larger
                maxTime = max(maxTime, ptimes{ch, bl}(end));
            end
        end
    end
    
    % Print the total time of the TCD file
    disp(['Total time of the TCD file: ', num2str(maxTime), ' seconds']);
