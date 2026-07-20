function varargout = process_psa_mx2( varargin )
% ASSUME: Inputs A and B have their time axes aligned (synchronised and bounds are the same).
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
sProcess.Comment     = 'Mean Flow Index NIRS';
sProcess.FileTag     = 'Mx';
sProcess.Category    = 'File2';
sProcess.SubGroup    = 'Frequency';
sProcess.Index       = 3; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'matrix', 'matrix'};
sProcess.nInputs     = 2;
sProcess.nMinFiles   = 1;
sProcess.nOutputs    = 1;

% Definition of the options
sProcess.options.abp_channel.Comment = 'ABP Channel: ';
sProcess.options.abp_channel.Type    = 'text';
sProcess.options.abp_channel.Value   = '';

sProcess.options.nirs_cardiac_event.Comment = 'NIRS cardiac event: ';
sProcess.options.nirs_cardiac_event.Type    = 'text';
sProcess.options.nirs_cardiac_event.Value   = 'cardiac';

sProcess.options.aux_cardiac_event.Comment = 'AUX cardiac event: ';
sProcess.options.aux_cardiac_event.Type    = 'text';
sProcess.options.aux_cardiac_event.Value   = 'cardiac';

sProcess.options.pct_filter.Comment = 'Percentage of extreme values to filter: ';
sProcess.options.pct_filter.Type    = 'value';
sProcess.options.pct_filter.Value   = {5, '%', 0};

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs_nirs, sInputs_aux) 
OutputFiles = {};

for iInput=1:length(sInputs_nirs)
    % Load NIRS recordings
    if strcmp(sInputs_nirs(iInput).FileType, 'data')     % Imported data structure
        sDataIn_nirs = in_bst_data(sInputs_nirs(iInput).FileName);
        events_nirs = sDataIn_nirs.Events;
    elseif strcmp(sInputs_nirs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn_nirs = in_bst(sInputs_nirs(iInput).FileName, [], 1, 1, 'no');
        sDataRaw_nirs = in_bst_data(sInputs_nirs(iInput).FileName, 'F');
        events_nirs = sDataRaw_nirs.F.events;
    end
    time_nirs = sDataIn_nirs.Time;
    channels_nirs = in_bst_channel(sInputs_nirs(iInput).ChannelFile);
    channel_mask_nirs = find_good_channel(sDataIn_nirs, channels_nirs, 'NIRS');

    % Load AUX recordings
    if strcmp(sInputs_aux(iInput).FileType, 'data')     % Imported data structure
        sDataIn_aux = in_bst_data(sInputs_aux(iInput).FileName);
        events_aux = sDataIn_aux.Events;
    elseif strcmp(sInputs_aux(iInput).FileType, 'raw')  % Continuous data file
        sDataIn_aux = in_bst(sInputs_aux(iInput).FileName, [], 1, 1, 'no');
        sDataRaw_aux = in_bst_data(sInputs_aux(iInput).FileName, 'F');
        events_aux = sDataRaw_aux.F.events;
    end

    % Check that time axes of AUX and NIRS are aligned
    if any(abs(sDataIn_aux.Time([1 end]) - sDataIn_nirs.Time([1 end])) > diff(sDataIn_nirs.Time(1:2)))
        error('Inconsistent time axes');
    end

    channel_bp = find_good_channel(sDataIn_aux, in_bst_channel(sInputs_aux(iInput).ChannelFile),...
                                   sProcess.options.abp_channel.Value);
    nb_bp = sum(channel_bp);
    if nb_bp == 0
        error('No BP channel available');
    elseif nb_bp > 1
        error('Multiple BP channels available');
    end
    abp_sig_idx_chan = find(channel_bp);
    
    bp_signal = sDataIn_aux.F(abp_sig_idx_chan, :);
    time_aux = sDataIn_aux.Time;
    sfreq_aux = 1/diff(time_aux(1:2));
    %TODO expose discard event label in options
    %TODO dedicated process that set discarded segments to nan
    %TODO decompose Mx: procude b2b AUX and NIRS signals with the same time
    %axis in the same data file with another channel definition
    % Then use brainstorm's correlation process as done in TCD
    idiscard_evt_aux = ismember({events_aux.label}, {'discard_BP', 'discard_cardiac'});
    time_mask_aux = true(1, size(time_aux, 2));
    if any(idiscard_evt_aux)      
       for evt_idx=find(idiscard_evt_aux)
           for isegment=1:size(events_aux(evt_idx).times, 2)
               event_times = events_aux(evt_idx).times(:, isegment);
               samplesBounds_aux = [round((event_times(1)-time_aux(1)) * sfreq_aux) + 1, ...
                                    round((event_times(2)-time_aux(1)) * sfreq_aux)];
               time_mask_aux(:, samplesBounds_aux(1):samplesBounds_aux(2)) = 0;
           end
       end
    end
    if sum(time_mask_aux) / length(time_mask_aux) < 0.5
        error('More than 50% of BP signal is discarded');
    end

    time_mask_nirs = discard_channels_time_mask(events_nirs, {channels_nirs.Channel.Name}, sDataIn_nirs.Time);
    time_mask_nirs(~channel_mask_nirs, :) = false;
    
    idiscard_evt_nirs = ismember({events_aux.label}, {'discard'});
    if any(idiscard_evt_nirs)      
       for evt_idx=find(idiscard_evt_nirs)
           for isegment=1:size(events_aux(evt_idx).times, 2)
               event_times = events_nirs(evt_idx).times(:, isegment);
               samplesBounds_nirs = [round((event_times(1)-time_nirs(1)) * sfreq_nirs) + 1, ...
                                    round((event_times(2)-time_nirs(1)) * sfreq_nirs)];
               time_mask_aux(:, samplesBounds_nirs(1):samplesBounds_nirs(2)) = 0;
           end
       end
    end

    icardiac_evt_aux = strcmp(sProcess.options.aux_cardiac_event.Value, {events_aux.label});
    icardiac_evt_nirs = strcmp(sProcess.options.nirs_cardiac_event.Value, {events_nirs.label});
    cardiac_events_aux = events_aux(icardiac_evt_aux).times;
    cardiac_events_nirs = events_nirs(icardiac_evt_nirs).times;
    if length(cardiac_events_aux) ~= length(cardiac_events_nirs)
        error('Inconsistent number of cardiac events between nirs (%d) and aux (%d)', ...
              length(cardiac_events_nirs), length(cardiac_events_aux));
    end
    idx_chans = find(channel_mask_nirs);
    mx_values = nan(size(channel_mask_nirs));
    for ichan=1:length(idx_chans)
        idx_chan = idx_chans(ichan);
        nirs_signal = sDataIn_nirs.F(idx_chan, :);
        chan_name = channels_nirs.Channel(idx_chan).Name;
        mx_out = Compute(bp_signal, time_aux, time_mask_aux, ...
                         nirs_signal, time_nirs, time_mask_nirs(idx_chan, :),...
                         events_nirs(icardiac_evt_nirs).times, ...
                         sProcess.options.pct_filter.Value{1}, ...
                         chan_name);
        mx_values(idx_chan) = mx_out;
    end   
    
    sDataOut = db_template('data');
    sDataOut.F            = mx_values;
    sDataOut.Comment      = [sInputs_nirs(iInput).Comment  ' | Mx' ];
    sDataOut.ChannelFlag  = sDataIn_nirs.ChannelFlag;
    sDataOut.Time         = [1];
    sDataOut.History      = sDataIn_nirs.History;
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.DisplayUnits = 'Pearson r';

    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs_nirs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_mx');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs_nirs(iInput).iStudy, OutputFile, sDataOut);

    OutputFiles{iInput} = OutputFile;
end

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
    if ~isempty(discard_channels_idx_events)
        assert(abs(masked_duration-nb_masked_samples/sfreq) < 1/sfreq);
        total_discard_duration = sum(diff(events(discard_channels_idx_events).times));
        total_masked_duration = sum(any(~time_mask)) / sfreq;
        assert(abs(total_discard_duration - total_masked_duration) < 1/sfreq);
    end
end

end


function label = protect_field_label(label)
label = strrep(label, '.', '_dot_');
end

function channel_mask = find_good_channel(sData, channel_mat, channel_filter)
channel_mask = false(size(sData.ChannelFlag));
nirs_idx_chan = channel_find(channel_mat.Channel, channel_filter);
channel_mask(nirs_idx_chan) = true;
channel_mask(sData.ChannelFlag ~= 1) = false;
end

%% ===== Compute =====
function [Mx] = Compute(bp_signal, bp_time, bp_time_mask, nirs_signal, nirs_time, nirs_time_mask, cardiac_events, pct_filter, channel_name)
% ASSUME bp_time and nirs_time are aligned

[b2b_nirs, nirs_beat_times] = apply_on_epochs(nirs_signal, nirs_time, nirs_time_mask, cardiac_events, @mean_amplitude);
[b2b_bp, bp_beat_times] = apply_on_epochs(bp_signal, bp_time, bp_time_mask, cardiac_events, @mean_amplitude);

assert(all(abs(bp_beat_times-nirs_beat_times) < 1e-6));

nan_mask = isnan(b2b_bp) | isnan(b2b_nirs);
b2b_nirs(nan_mask) = [];
b2b_bp(nan_mask) = [];

[b2b_bp_filtered, b2b_nirs_filtered] = pct_filter2(b2b_bp, b2b_nirs, pct_filter);

[r, p] = corrcoef(b2b_bp_filtered, b2b_nirs_filtered);

if p > 0.05
    warning('%s: No significant correlation between block signals of BP and CBFi', channel_name);
end

Mx = r(1,2);
end

function event_mask = raster_events(events, time)
event_mask = false(size(time));
if ~isempty(events)
    samples = time_to_sample_idx(events, time);
    event_mask(samples) = 1;
end
end

function [values, times] = apply_on_epochs(signal, signal_time, time_mask, epoch_times, func, channel_name)

sampling_freq = 1 ./ diff(signal_time(1:2));
epoch_samples = round((epoch_times-signal_time(1)) .* sampling_freq) + 1;

values = NaN(size(signal,1), length(epoch_samples)-1);
times = zeros(1, length(epoch_samples)-1);
for iepoch=1:(length(epoch_samples)-1)
    [value, time] = func(signal(:, epoch_samples(iepoch):epoch_samples(iepoch+1)), ...
                         signal_time(epoch_samples(iepoch):epoch_samples(iepoch+1)));
    if all(time_mask(epoch_samples(iepoch):epoch_samples(iepoch+1)))
        values(:, iepoch) = value;
    % else
    %     fprintf('Channel %s: Ignore cardiac epoch from %1.1fs to %1.1fs because it contains a discard segment\n', ...
    %             channel_name, signal_time(epoch_samples(iepoch)), signal_time(epoch_samples(iepoch+1)));
    end
    times(1, iepoch) = time;
end

end

function [mean_value, time] = mean_amplitude(heart_beat_signal, heart_beat_time)
    mean_value = mean(heart_beat_signal,2);
    time = heart_beat_time(1) + (heart_beat_time(end) - heart_beat_time(1)) / 2;
end

function [sig1_filtered, sig2_filtered] = pct_filter2(sig1, sig2, pct)
pct = pct / 100;
[~, sig1_sorted_idx] = sort(sig1);
sig1_low_idx = max(1, round(length(sig1)*pct/2));
sig1_high_idx = length(sig1) - sig1_low_idx;

sig1 = sig1(sig1_sorted_idx(sig1_low_idx:sig1_high_idx));
sig2 = sig2(sig1_sorted_idx(sig1_low_idx:sig1_high_idx));

[~, sig2_sorted_idx] = sort(sig2);
sig2_low_idx = max(1, round(length(sig2)*pct/2));
sig2_high_idx = length(sig2) - sig2_low_idx;

sig1_filtered = sig1(sig2_sorted_idx(sig2_low_idx:sig2_high_idx));
sig2_filtered = sig2(sig2_sorted_idx(sig2_low_idx:sig2_high_idx));

end


function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end
