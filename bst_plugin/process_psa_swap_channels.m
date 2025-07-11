function varargout = process_psa_swap_channels( varargin )

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
% Authors: Thomas Vincent (2025-)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Swap channels';
sProcess.FileTag     = 'cswap';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Standardize';
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

sProcess.options.channelname_1.Comment = 'Name of first channel to swap:';
sProcess.options.channelname_1.Type    = 'text';
sProcess.options.channelname_1.Value   = '';

sProcess.options.channelname_2.Comment = 'Name of second channel to swap:';
sProcess.options.channelname_2.Type    = 'text';
sProcess.options.channelname_2.Value   = '';
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
    idx_chan1 = channel_find(channels.Channel, sProcess.options.channelname_1.Value);
    idx_chan2 = channel_find(channels.Channel, sProcess.options.channelname_2.Value);

    swapped_chan_order = 1:length(channels.Channel);
    swapped_chan_order([idx_chan1 idx_chan2]) = swapped_chan_order([idx_chan2 idx_chan1]);

    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = sDataIn.F(swapped_chan_order, :);
    sDataOut.Comment      = [sDataIn.Comment ' | Swapped Chans'];
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time;
    sDataOut.History      = sDataIn.History;
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.Events       = events;
    sDataOut.DisplayUnits = sDataIn.DisplayUnits;
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_cswap');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    OutputFiles{iInput} = OutputFile;
end
end


function [values, times] = apply_on_epochs(signal, signal_time, epoch_times, func)

sampling_freq = 1 ./ diff(signal_time(1:2));
epoch_samples = round(epoch_times .* sampling_freq);
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
    [systolic_peak, i_sys_peak] = max(heart_beat_signal);
    diastolic_peak = min(heart_beat_signal(i_sys_peak:end, :));
    pulsatility_value = (systolic_peak - diastolic_peak) ./ diastolic_peak;
    time = heart_beat_time(1) + (heart_beat_time(end) - heart_beat_time(1)) / 2;
end