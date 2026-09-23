function varargout = process_psa_curate_event( varargin )

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
% Authors: Thomas Vincent (2026-)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
%TOCHECK: how do we limit the input file types (only NIRS data)?
sProcess.Comment     = 'Curate event';
sProcess.FileTag     = 'curated';
sProcess.Category    = 'Filter';
sProcess.SubGroup    = 'Pre-process';
sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in
%                             the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};

% Definition of the options
sProcess.options.channelnames.Comment = 'Channel names (Coma-separated values):';
sProcess.options.channelnames.Type    = 'text';
sProcess.options.channelnames.Value   = '';

sProcess.options.discard_global_events.Comment = 'Discard events (apply to all channels, coma-separated):';
sProcess.options.discard_global_events.Type    = 'text';
sProcess.options.discard_global_events.Value   = ''; 

sProcess.options.use_discard_channel_events.Comment = 'Use discard_<channel> events to discard channel-specific segments';
sProcess.options.use_discard_channel_events.Type       = 'checkbox';
sProcess.options.use_discard_channel_events.Value      = 1;

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function sInputs = Run(sProcess, sInputs) %#ok<DEFNU>
if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs.FileName, 'Events');
    events = sDataIn.Events;
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs.FileName, [], 1, 1, 'no');
    sDataRaw = in_bst_data(sInputs.FileName, 'F');
    events = sDataRaw.F.events;
else
    error('Unhandled file type %s', sInputs.FileType);
end
if sProcess.options.use_discard_channel_events.Value
    channel_discard_event_prefixes = {'discard_'};
else
    channel_discard_event_prefixes = {};
end
chanMat = in_bst_channel(sInputs.ChannelFile);
sInputs.A = Compute(sInputs.A, sInputs.TimeVector, events, ...
                    chanMat.Channel, ...
                    sProcess.options.channelnames.Value, ...
                    sProcess.options.discard_global_events.Value, ...
                    channel_discard_event_prefixes);
end


function signal_curated = Compute(signal, time, events, channels, channel_filter, global_events, channel_discard_event_prefixes) 

channel_names = {channels.Name};

nb_channels = size(channels, 2);
if ~isempty(channel_filter)
    idx_chans = channel_find(channels, channel_filter);
    chan_mask = false(1, nb_channels);
    chan_mask(idx_chans) = 1;
else
    chan_mask = true(1, nb_channels);
end

    
% Filter channel-specific discard segments
if ~isempty(channel_discard_event_prefixes)
    for ievent=1:length(channel_discard_event_prefixes)
        event_prefix = channel_discard_event_prefixes{ievent};
        if ievent==1
            time_mask = channels_time_mask(events, channel_names, time, event_prefix);
        else
            time_mask = time_mask & channels_time_mask(events, channel_names, time, event_prefix);
        end
    end
else
    time_mask = true(nb_channels, size(time, 2));
end

% Filter global discard segments
if ~isempty(global_events)
   discard_event_names = cellfun(@strtrim, strsplit(global_events, ','),...
                                 'UniformOutput', 0);
   sfreq = 1/diff(time(1:2));
   global_time_mask = true(size(time, 2), 1);
   for idiscard_event=1:length(discard_event_names)
        idiscard_evt = find(strcmp(discard_event_names{idiscard_event}, {events.label}));
        if ~isempty(idiscard_evt)
            event_times = events(idiscard_evt).times;
            for evt_idx=1:size(event_times, 2)
               for isegment=1:size(event_times, 2)
                   time_bounds = event_times(:, isegment);
                   samplesBounds = [round((time_bounds(1)-time(1)) * sfreq) + 1, ...
                                    round((time_bounds(2)-time(1)) * sfreq)];
                   global_time_mask(samplesBounds(1):samplesBounds(2)) = false;
               end
            end
        end
   end
   % prop_discard = 1 - sum(global_time_mask)/size(sDataIn.Time, 2);
   % fprintf('Percent of signal discarded (non-channel-specific) : %1.0f%%\n', ...
   %         prop_discard * 100);

   for ichan=1:nb_channels
        time_mask(ichan, :) = time_mask(ichan, :) & global_time_mask';
   end
end
time_mask(~chan_mask, :) = false;

signal_curated = signal;
signal_curated(~time_mask) = nan;

end

function time_mask=channels_time_mask(events, channel_names, time, event_prefix)
pattern = [event_prefix '.*'];
time_mask = true(length(channel_names), size(time, 2));
discard_idx_events = find(~cellfun(@isempty, regexpi({events.label}, pattern)));
if ~isempty(discard_idx_events)
    channel_map = containers.Map(channel_names, 1:length(channel_names));
    sfreq = 1/diff(time(1:2));
    nb_masked_samples = 0;
    masked_duration = 0;
    discard_channels_idx_events = [];
    for evt_idx=discard_idx_events
        suffix = events(evt_idx).label((length(event_prefix)+1):end);
        if channel_map.isKey(suffix)
            ichannel = channel_map(suffix);
            for isegment=1:size(events(evt_idx).times, 2)
                event_times = events(evt_idx).times(:, isegment);
                samplesBounds = [round((event_times(1)-time(1)) * sfreq) + 1, ...
                                 round((event_times(2)-time(1)) * sfreq)];
                time_mask(ichannel, samplesBounds(1):samplesBounds(2)) = false;
                nb_masked_samples = nb_masked_samples + diff(samplesBounds) + 1;
                masked_duration = masked_duration + diff(events(evt_idx).times(:, isegment));
                if 0
                    fprintf('Discard segment [%1.2f-%1.2f]=%1.2f s, nb samples=%d -> %1.2f s\n', ...
                            events(evt_idx).times(1, isegment), events(evt_idx).times(1, isegment),...
                            masked_duration, nb_masked_samples, nb_masked_samples / sfreq);
                end
            end
            discard_channels_idx_events = [discard_channels_idx_events evt_idx];
        end
    end
    if 0 && ~isempty(discard_channels_idx_events)
        assert(abs(masked_duration-nb_masked_samples/sfreq) < 1/sfreq);
        total_discard_duration = sum(diff(events(discard_channels_idx_events).times));
        total_masked_duration = sum(any(~time_mask)) / sfreq;
        assert(abs(total_discard_duration - total_masked_duration) < 1/sfreq);
    end
end

end

