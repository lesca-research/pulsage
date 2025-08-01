function varargout = process_psa_average_time_robust( varargin )

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
% Authors: Aymen Zire, Thomas Vincent (2025-)
%
% Adapted from process_average_time.m

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Average time (robust)';
sProcess.FileTag     = @GetFileTag;
sProcess.Category    = 'Filter';
sProcess.SubGroup    = 'Average';
sProcess.Index       = 304; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options

sProcess.options.pct_filter_range.Comment = 'Percentile filter range (percent):';
sProcess.options.pct_filter_range.Type    = 'range';
sProcess.options.pct_filter_range.Value   = {[2.5 97.5], '%', 1};

sProcess.options.discard_bad_events.Comment = 'Use bad_ events to discard segments';
sProcess.options.discard_bad_events.Type       = 'checkbox';
sProcess.options.discard_bad_events.Value      = 0;

% === FUNCTION
sProcess.options.label2.Comment = '<U><B>Function</B></U>:';
sProcess.options.label2.Type    = 'label';
sProcess.options.avg_func.Comment = {'Arithmetic average:  <FONT color="#777777">mean(x)</FONT>', ...
                                     'Root mean square (RMS):  <FONT color="#777777">sqrt(sum(x.^2)/N)</FONT>', ...
                                     'Standard deviation:  <FONT color="#777777">sqrt(var(x))</FONT>', ...
                                     'Median:  <FONT color="#777777">median(x)</FONT>', ...
                                     'Minimum:  <FONT color="#777777">min(x)</FONT>', ...
                                     'Max:  <FONT color="#777777">max(x)</FONT>'; ...
                                     'mean', 'rms', 'std', 'median', 'min', 'max'};
sProcess.options.avg_func.Type    = 'radio_label';
sProcess.options.avg_func.Value   = 'mean';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
    % Old version of the process: option isstd={0,1}
    if isfield(sProcess.options, 'isstd') && ~isempty(sProcess.options.isstd) && sProcess.options.isstd.Value
        fileTag = 'std';
    % New version of the process: avg_fun={'mean','std','rms','median'}
    elseif isfield(sProcess.options, 'avg_func') && ~isempty(sProcess.options.avg_func) && ~isempty(sProcess.options.avg_func.Value)
        fileTag = sProcess.options.avg_func.Value;
    else
        fileTag = 'mean';
    end
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
channels = in_bst_channel(sInput.ChannelFile);
channel_names = {channels.Channel.Name};
nb_channels = size(sInput.A, 1);
time_mask = true(nb_channels, size(sInput.TimeVector, 2));

if ~ismatrix(sInput.A)
    error('Unsupported number of dimensions');
end

% Filter bad segments
if sProcess.options.discard_bad_events.Value
    sEvents = in_bst_data(sInput.FileName, 'Events');
    bad_idx_events = find(~cellfun(@isempty, regexpi({sEvents.Events.label}, ['bad_.*'])));
    if ~isempty(bad_idx_events)
        channel_map = containers.Map(channel_names, 1:length(channel_names));
        sfreq = 1/diff(sInput.TimeVector(1:2));
        for evt_idx=bad_idx_events
            suffix = sEvents.Events(evt_idx).label(5:end);
            if channel_map.isKey(suffix)
                ichannel = channel_map(suffix);
                samplesBounds = round((sEvents.Events(evt_idx).times-sInput.TimeVector(1)) * sfreq);
                time_mask(ichannel, samplesBounds(1):samplesBounds(2)) = 0;
                fprintf('Discard channel-specific bad segment(s) using event %s\n', sEvents.Events(evt_idx).label);
            else
                fprintf('Discard bad segment(s) for all channels using event %s\n', sEvents.Events(evt_idx).label);
                time_mask(:, samplesBounds(1):samplesBounds(2)) = 0;
            end
        end
    end
end

A = zeros(size(sInput.A, 1), 1);
operation = GetFileTag(sProcess);
for ichan=1:nb_channels
    channel_data = sInput.A(ichan, time_mask(ichan, :), :);
    time = sInput.TimeVector(time_mask(ichan, :));
    nb_samples = length(channel_data);
    low = round(sProcess.options.pct_filter_range.Value{1}(1)/100 * nb_samples);
    high = round(sProcess.options.pct_filter_range.Value{1}(2)/100 * nb_samples);
    [sorted_data, sort_idx] = sort(channel_data);
    time = time(sort_idx);
    channel_data = sorted_data(low:high);
    time = time(low:high);
    if any(isnan(channel_data))
        warning('Nan value(s) in filtered data of channel "%s"\n', channel_names{ichan})
    end

    % Apply function
    result_itime = 1;
    switch operation
        case 'mean'
            result = mean(channel_data, 2);
        case 'rms'
            result = sqrt(sum(channel_data.^2, 2) / length(channel_data));
        case 'std'
            result = sqrt(var(channel_data, 0, 2));
        case 'median'
            result = median(channel_data, 2);
        case 'min'
            [result, result_itime] = min(channel_data, 2);
        case 'max'
            [result, result_itime] = max(channel_data, 2);
    end
    if any(isnan(result))
        warning('Nan value(s) in computation of %s for channel "%s"\n', operation, channel_names{ichan})
    end
    A(ichan, :, :) = result;
end
% Copy values to represent the time window
sInput.A = [A, A];
% Keep only first and last time values

if (length(sInput.TimeVector) > 2)
    sInput.TimeVector = time(result_itime) + [0, sInput.TimeVector(2)-sInput.TimeVector(1)];
else
    sInput.TimeVector = time(result_itime) + [0, 1e-6];
end

% Build file tag
sInput.CommentTag = [GetFileTag(sProcess) '(' process_extract_time('GetTimeString',sProcess,sInput) ')'];
% Do not keep the Std/TFmask fields in the output
if isfield(sInput, 'Std') && ~isempty(sInput.Std)
    sInput.Std = [];
end
if isfield(sInput, 'TFmask') && ~isempty(sInput.TFmask)
    sInput.TFmask = [];
end

end



