function varargout = process_psa_pulsatility( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Thomas Vincent (2024-)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
%TOCHECK: how do we limit the input file types (only NIRS data)?
sProcess.Comment     = 'Pulsatility';
sProcess.FileTag     = 'pulsatility';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Pre-process';
sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in
%                             the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options

sProcess.options.channelnames.Comment = 'Channel names (Coma-separated values):';
sProcess.options.channelnames.Type    = 'text';
sProcess.options.channelnames.Value   = '';

sProcess.options.heart_beat_event.Comment = 'Heart beat event';
sProcess.options.heart_beat_event.Type    = 'text';
sProcess.options.heart_beat_event.Value   = 'cardiac';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};

for iInput=1:length(sInputs)
    % Load recordings
    if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iInput).FileName);
        events = sDataIn.Events;
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iInput).FileName, 'F');
        events = sDataRaw.F.events;
    end
    channels = in_bst_channel(sInputs(iInput).ChannelFile);
    nb_channels = size(channels.Channel, 2);
    if ~isempty(sProcess.options.channelnames.Value)
        idx_chans = channel_find(channels.Channel, sProcess.options.channelnames.Value);
        chan_mask = false(1, nb_channels);
        chan_mask(idx_chans) = 1;
    else
        chan_mask = true(1, nb_channels);
    end

    signal = sDataIn.F(chan_mask, :)';
    
    ievent = strcmp(sProcess.options.heart_beat_event.Value, ...
                    {events.label});

    nb_selected_chans = sum(chan_mask);
    pulsatility = zeros(length(events(ievent).times(1,:)) - 1, nb_selected_chans);
    chan_indices = find(chan_mask);
    for iChan = 1:nb_selected_chans
        sig_chan = sDataIn.F(chan_indices(iChan), :)';  % signal 1D pour le canal i
        [puls_chan, pulsatility_time] = apply_on_epochs(sig_chan, sDataIn.Time,...
                                                        events(ievent).times(1,:),...
                                                        @pulsatility);
        pulsatility(:, iChan) = puls_chan;
    end


    iSamples = round(pulsatility_time / diff(sDataIn.Time(1:2)));
    signal_full = sDataIn.F(:, iSamples);
    signal_full(chan_mask, :) = pulsatility';

    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = signal_full;
    sDataOut.Comment      = [sDataIn.Comment ' | Pulsatility'];
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = pulsatility_time;
    sDataOut.History      = sDataIn.History;
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.Events       = events;
    sDataOut.DisplayUnits = sDataIn.DisplayUnits;
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_mavg');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    OutputFiles{iInput} = OutputFile;
end
end


function [values, times] = apply_on_epochs(signal, signal_time, epoch_times, func)

sampling_freq = 1 ./ diff(signal_time(1:2));
epoch_samples = round((epoch_times-signal_time(1)) .* sampling_freq) + 1;

values = zeros(length(epoch_samples)-1, size(signal,2));
times = zeros(1, length(epoch_samples)-1);
for iepoch=1:(length(epoch_samples)-1)
    [value, time] = func(signal(epoch_samples(iepoch):epoch_samples(iepoch+1), :), ...
                         signal_time(epoch_samples(iepoch):epoch_samples(iepoch+1)));
    values(iepoch,:) = value;
    times(1, iepoch) = time;
end

end

function [pulsatility_value, time] = pulsatility(heart_beat_signal, heart_beat_time)
    systolic_peak = max(heart_beat_signal);
    diastolic_peak = min(heart_beat_signal);
    pulsatility_value = abs(systolic_peak - diastolic_peak) ./ systolic_peak;
    time = heart_beat_time(1) + (heart_beat_time(end) - heart_beat_time(1)) / 2;
end
