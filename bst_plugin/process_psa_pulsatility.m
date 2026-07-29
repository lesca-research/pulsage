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
sProcess.SubGroup    = 'PulsAge';
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

sProcess.options.flip_sign.Comment = 'Flip sign';
sProcess.options.flip_sign.Type    = 'checkbox';
sProcess.options.flip_sign.Value   = 0;


sProcess.options.label2.Comment = '<U><B>Function to apply to each epoch signal</B></U>:';
sProcess.options.label2.Type    = 'label';
sProcess.options.function.Comment = {'Pulse CBF', ...
                                     'Time to Pulse CBF', ...
                                     'Pulse Amplitude:  <FONT color="#777777">(systole - diastole)/diastole</FONT>', ...
                                     'Heart Rate:  <FONT color="#777777">1/duration*60</FONT>'; ...
                                     'pCBF', 'tt_pCBF', 'pulse_amplitude', 'heart_rate'};
sProcess.options.function.Type    = 'radio_label';
sProcess.options.function.Value   = 'pulse_amplitude_normed';

sProcess.options.use_discard_channel_events.Comment = 'Use discard_<channel> events to discard channel-specific segments';
sProcess.options.use_discard_channel_events.Type       = 'checkbox';
sProcess.options.use_discard_channel_events.Value      = 1;

sProcess.options.discard_global_events.Comment = 'Discard events (apply to all channels, coma-separated):';
sProcess.options.discard_global_events.Type    = 'text';
sProcess.options.discard_global_events.Value   = ''; 

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
    channel_names = {channels.Channel.Name};

    nb_channels = size(channels.Channel, 2);
    if ~isempty(sProcess.options.channelnames.Value)
        idx_chans = channel_find(channels.Channel, sProcess.options.channelnames.Value);
        chan_mask = false(1, nb_channels);
        chan_mask(idx_chans) = 1;
    else
        chan_mask = true(1, nb_channels);
    end

    
    % Filter channel-specific discard segments
    if sProcess.options.use_discard_channel_events.Value
        time_mask = discard_channels_time_mask(events, channel_names, sDataIn.Time);
    else
        time_mask = true(nb_channels, size(sDataIn.Time, 2));
    end

    % Filter global discard segments
    if ~isempty(sProcess.options.discard_global_events.Value)
       discard_event_names = cellfun(@strtrim, strsplit(sProcess.options.discard_global_events.Value, ','),...
                                     'UniformOutput', 0);
       sfreq = 1/diff(sDataIn.Time(1:2));
       global_time_mask = true(size(sDataIn.Time, 2), 1);
       for idiscard_event=1:length(discard_event_names)
            idiscard_evt = find(strcmp(discard_event_names{idiscard_event}, {events.label}));
            if ~isempty(idiscard_evt)
                event_times = events(idiscard_evt).times;
                for evt_idx=1:size(event_times, 2)
                   for isegment=1:size(event_times, 2)
                       time_bounds = event_times(:, isegment);
                       samplesBounds = [round((time_bounds(1)-sDataIn.Time(1)) * sfreq) + 1, ...
                                        round((time_bounds(2)-sDataIn.Time(1)) * sfreq)];
                       global_time_mask(samplesBounds(1):samplesBounds(2)) = false;
                   end
                end
            end
       end
       prop_discard = 1 - sum(global_time_mask)/size(sDataIn.Time, 2);
       fprintf('Percent of signal discarded (non-channel-specific) : %1.0f%%\n', ...
               prop_discard * 100);
    
       for ichan=1:nb_channels
            time_mask(ichan, :) = time_mask(ichan, :) & global_time_mask';
       end
    end
    time_mask(~chan_mask, :) = false;

    ievent = strcmp(sProcess.options.heart_beat_event.Value, ...
                    {events.label});

    pulsatility = nan(length(events(ievent).times(1,:)) - 1, nb_channels);
    for chan_idx = 1:nb_channels
        if sum(time_mask(chan_idx, :)) / size(sDataIn.Time, 2) > 0.5
            sig_chan = sDataIn.F(chan_idx, :)';
            time_mask_chan = time_mask(chan_idx, :);
            if sProcess.options.flip_sign.Value
                sig_chan = sig_chan * -1;
            end
    
            switch sProcess.options.function.Value
                case 'pCBF'
                    func = @gradient_max;
                case 'tt_pCBF'
                    func = @gradient_argmax;
                case 'pulse_amplitude'
                    func = @pulse_amplitude;
                case 'pulse_amplitude_normed'
                    func = @pulse_amplitude_normed;
                case 'heart_rate'
                    func = @heart_rate;
            end
            [puls_chan, pulsatility_time] = apply_on_epochs(sig_chan, sDataIn.Time, ...
                                                            time_mask_chan, ...
                                                            events(ievent).times(1,:),...
                                                            func, channel_names{chan_idx});
            pulsatility(:, chan_idx) = puls_chan;
        end
    end

    %iSamples = round((pulsatility_time - sDataIn.Time(1))  / diff(sDataIn.Time(1:2)));
    %signal_full = sDataIn.F(:, iSamples);
    %signal_full(chan_mask, :) = pulsatility';

    % Filter extened events to align with new time:
    for ievt_group=1:length(events)
        if size(events(ievt_group).times, 1) == 2
            m_out_left = events(ievt_group).times(1, :) < pulsatility_time(1);
            events(ievt_group).times(1, m_out_left) = pulsatility_time(1);

            m_out_right = events(ievt_group).times(2, :) > pulsatility_time(end);
            events(ievt_group).times(1, m_out_right) = pulsatility_time(end);
    
            m_empty = diff(events(ievt_group).times) == 0;
            if any(m_empty)
                events(ievt_group).times(:, m_empty) = [];
                events(ievt_group).epochs(m_empty) = [];
            end
        end
    end
    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = pulsatility';
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


function [values, times] = apply_on_epochs(signal, signal_time, time_mask, epoch_times, func, channel_name)

sampling_freq = 1 ./ diff(signal_time(1:2));
epoch_samples = round((epoch_times-signal_time(1)) .* sampling_freq) + 1;

values = NaN(length(epoch_samples)-1, size(signal,2));
times = zeros(1, length(epoch_samples)-1);
for iepoch=1:(length(epoch_samples)-1)
    if epoch_samples(iepoch) == epoch_samples(iepoch+1)
        error('Epoch too short at %f s', epoch_samples(iepoch)/sampling_freq + signal_time(1))
    end
    [value, time] = func(signal(epoch_samples(iepoch):epoch_samples(iepoch+1), :), ...
                         signal_time(epoch_samples(iepoch):epoch_samples(iepoch+1)));
    if all(time_mask(epoch_samples(iepoch):epoch_samples(iepoch+1)))
        values(iepoch,:) = value;
    % else
    %     fprintf('Channel %s: Ignore cardiac epoch from %1.1fs to %1.1fs because it contains a discard segment\n', ...
    %             channel_name, signal_time(epoch_samples(iepoch)), signal_time(epoch_samples(iepoch+1)));
    end
    times(1, iepoch) = time;
end

end

function [pulsatility_value, time] = pulse_amplitude_normed(heart_beat_signal, heart_beat_time)
    systolic_peak = max(heart_beat_signal);
    diastolic_peak = min(heart_beat_signal);

    % resistance index: systolic_peak - diastolic_peak ./ systolic_peak
    pulsatility_value = (systolic_peak - diastolic_peak) ./ mean(heart_beat_signal);
    time = heart_beat_time(1) + (heart_beat_time(end) - heart_beat_time(1)) / 2;
end

function [pulsatility_value, time] = pulse_amplitude(heart_beat_signal, heart_beat_time)
    systolic_peak = max(heart_beat_signal);
    diastolic_peak = min(heart_beat_signal);

    % resistance index: systolic_peak - diastolic_peak ./ systolic_peak
    % pulsatility_value = (systolic_peak - diastolic_peak) ./ mean(heart_beat_signal);
    pulsatility_value = (systolic_peak - diastolic_peak);
    time = heart_beat_time(1) + (heart_beat_time(end) - heart_beat_time(1)) / 2;
end

function [heart_rate, time] = heart_rate(heart_beat_signal, heart_beat_time)
    heart_rate = 1/(heart_beat_time(end) - heart_beat_time(1))*60;
    time = heart_beat_time(1) + (heart_beat_time(end) - heart_beat_time(1)) / 2;
end

function [pulsatility_value, time] = gradient_max(heart_beat_signal, heart_beat_time)
    gradient = diff(heart_beat_signal);
    [pulsatility_value, max_itime] = max(gradient);
    time = heart_beat_time(max_itime);
end

function [pulsatility_value, time] = gradient_argmax(heart_beat_signal, heart_beat_time)
    gradient = diff(heart_beat_signal);
    [~, max_itime] = max(gradient);
    time = heart_beat_time(max_itime);
    pulsatility_value = time;
end


function time_mask=discard_channels_time_mask(events, channel_names, time)
time_mask = true(length(channel_names), size(time, 2));
discard_idx_events = find(~cellfun(@isempty, regexpi({events.label}, 'discard_.*')));
if ~isempty(discard_idx_events)
    channel_map = containers.Map(channel_names, 1:length(channel_names));
    sfreq = 1/diff(time(1:2));
    nb_masked_samples = 0;
    masked_duration = 0;
    discard_channels_idx_events = [];
    for evt_idx=discard_idx_events
        suffix = events(evt_idx).label(9:end);
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

