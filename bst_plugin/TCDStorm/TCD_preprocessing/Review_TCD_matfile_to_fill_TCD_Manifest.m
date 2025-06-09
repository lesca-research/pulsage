close all
clear all
clc

%% Code is written by Hanieh Mohammadi Version 19 Feb 2024
%% For any question please contact hanieh.bme@gmail.com or hanieh.mohammadi@polymtl.ca

%% Example of input
% pathname='/Users/haniehm/Documents/ACTIONcrNov1st2023/ACTION_TCD/';
% fname='ACR-0001-00230_TCD_T12_rest_staning_LMCA_RMCA.mat';

  fname = 'RC104_T0_rest_stand_TCD_RMCA_LMCA';
  pathname= '/Users/haniehm/Documents/ACTIONcrNov1st2023/ACTION_TCD/TCDStorm/TroubleShooting_Recardio/'
  %pathname = '/Users/haniehm/Documents/ACTIONcrNov1st2023/ACTION_TCD/';
  %csv_folder_path = '/Users/haniehm/Documents/ACTIONcrNov1st2023/ACTION_TCD/TCDStorm/TroubleShooting_Recardio/';
 % csv_filename = fullfile(csv_folder_path, 'TCD_manifest.csv');

load([pathname fname]);

if ~exist('data','var'),
    error('No data! Select a mat file that contains data and was created with Export Matlab or later (LabChart for Windows 7.2 or later)');
    return
end

[numchannels, numblocks] = size(datastart);

if numblocks == 2
    datastart = datastart(:, 2);
    dataend = dataend(:, 2);
    numblocks = 1;
end

% Modified plot creation logic
totalPlots = numchannels * numblocks;
plotCounter = 0;
ax = [];  % Array to store axes handles

for ch = 1:numchannels
    for bl = 1:numblocks
        if (datastart(ch,bl) ~= -1) % Empty blocks excluded
            plotCounter = plotCounter + 1;

            if mod(plotCounter, 2) == 1
                figureHandle = figure;
                scrsz = get(0,'ScreenSize');
                set(figureHandle, 'Position', scrsz([3 4 3 4]).*[1 1 6 6]/8);
                ax = [];  % Reset axes handle array for new figure
            end

            % Plotting
            pdata = data(datastart(ch,bl):dataend(ch,bl));
            if exist('scaleunits','var') % 16-bit data
                pdata = (pdata + scaleoffset(ch,bl)).* scaleunits(ch,bl);
            end
            ptime = [0:length(pdata)-1]/samplerate(ch,bl);

            % Adjusted subplot indexing
            ax(end+1) = subplot(2, 1, 2 - mod(plotCounter, 2));
            plot(ptime, pdata), hold on

            % Titles
            if bl == 1
                title(titles(ch,:));
            end

            % Ranges and units
            xlabel('Time (s)');
            xlim([0 max(ptime)]);
            ylabel(unittext(unittextmap(ch,bl),:));
            ylim([min(pdata) max(pdata)]);

            % Comments (if any)
            if exist('com','var'),
                temp = find((com(:,2) == bl & (com(:,1) == ch | com(:,1) == -1)));
                comtickpos = com(:,3);
                comtextmap = com(:,5);
                for m = 1:length(temp)
                    temp2 = temp(m);
                    x = round(comtickpos(temp2) * ...        % Convert from tick rate position
                              (samplerate(ch,bl) / tickrate(bl)));  % to sample rate position
                    t = x / samplerate(ch,bl);
                    plot([t t],[min(pdata) max(pdata)],'r:');
                    escapedComtext = strrep(comtext(comtextmap(temp2),:), '_', '\_');
                    text(t, min(pdata) + 0.1*(max(pdata) - min(pdata)), escapedComtext, ...
                        'Rotation', 90, 'BackgroundColor', 'white');
                end
            end

            % Link axes for simultaneous zooming and panning on x-axis
            if mod(plotCounter, 2) == 0 || plotCounter == totalPlots
                linkaxes(ax, 'x');
            end
        end
    end
end

%% Print comtext with positions
if exist('com','var') && exist('samplerate','var')
    fprintf('\nComments and their timestamps:\n');
    for i = 1:size(com, 1)
        tickpos = com(i, 3);
        comIndex = com(i, 5);
        timeInSeconds = tickpos / samplerate(1, 1); % Adjust if samplerate varies

        if comIndex > 0 && comIndex <= size(comtext, 1)
            currentComtext = strtrim(comtext(comIndex, :));
            if isempty(currentComtext) || all(isspace(currentComtext))
                fprintf('User did not input a trigger (such as rest start, rest end, standing start, or standing end) during data collection at %.3f s\n', timeInSeconds);
            else
                escapedComtext = strrep(currentComtext, '_', '_');
                fprintf('%s started at %.3f s\n', escapedComtext, timeInSeconds);
            end
        else
            fprintf('User did not input a trigger (such as rest start, rest end, standing start, or standing end) during data collection at %.3f s\n', timeInSeconds);
            % Print signal termination time
           

        end
    end
end

signalEndTime = max(ptime);
fprintf('Signal ends at %.3f s\n', signalEndTime);

% print signal end