function varargout = process_psa_spline_fit( varargin )

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
sProcess.Comment     = 'Spline fit';
sProcess.FileTag     = 'spline';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Pre-process';
sProcess.Index       = 1306; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in
%                             the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options

sProcess.options.channel_names.Comment = 'Channels (comma-separated): ';
sProcess.options.channel_names.Type    = 'text';
sProcess.options.channel_names.Value   = '';

sProcess.options.spline_channel_names.Comment = 'Output Channels (comma-separated, leave blank for inplace processing): ';
sProcess.options.spline_channel_names.Type    = 'text';
sProcess.options.spline_channel_names.Value   = '';


sProcess.options.ignore_gap_events.Comment = 'Ignore gap events (gap_<channel_name})';
sProcess.options.ignore_gap_events.Type    = 'checkbox';
sProcess.options.ignore_gap_events.Value   = 1;

sProcess.options.ignore_discard_events.Comment = 'Ignore discard events (discard_<channel_name})';
sProcess.options.ignore_discard_events.Type    = 'checkbox';
sProcess.options.ignore_discard_events.Value   = 1;

sProcess.options.prior_interpolation_dt.Comment = 'Dt for prior interpolation';
sProcess.options.prior_interpolation_dt.Type    = 'value';
sProcess.options.prior_interpolation_dt.Value   = {0.05, 'sec', 2};

sProcess.options.spline_knot_dt.Comment = 'Dt for knot temporal spacing';
sProcess.options.spline_knot_dt.Type    = 'value';
sProcess.options.spline_knot_dt.Value   = {0.15, 'sec', 2};
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};

if isempty(sProcess.options.spline_channel_names) && ~isempty(sProcess.options.channel_names)
    error('Channel names cannot be empty if destination spline channel names are provided');
end

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
    chanMat = in_bst_channel(sInputs(iInput).ChannelFile);
    nb_channels = size(chanMat.Channel, 2);
    chan_mask = (sDataIn.ChannelFlag == 1)'; 
    if ~isempty(sProcess.options.channel_names.Value)
        idx_chans = channel_find(chanMat.Channel, sProcess.options.channel_names.Value);
        chan_mask_tmp = false(1, nb_channels);
        chan_mask_tmp(idx_chans) = true;
        chan_mask = chan_mask & chan_mask_tmp;
    end
    channels = chanMat.Channel(chan_mask);

    signal = sDataIn.F(chan_mask, :)';
    nb_proc_channels = size(signal, 2);

    % Downsample / moving average
    interp_dt = sProcess.options.prior_interpolation_dt.Value{1};
    assert(interp_dt > diff(sDataIn.Time(1:2)));
    interp_time = sDataIn.Time(1):interp_dt:sDataIn.Time(end);
    if interp_time(end) ~= sDataIn.Time(end)
        interp_time = [interp_time sDataIn.Time(end)];
    end
    signal_interp = zeros(length(interp_time), nb_proc_channels);
    for ichan=1:nb_proc_channels
        signal_interp(:, ichan) = interp1(sDataIn.Time, signal(:, ichan)', interp_time, 'spline');
    end
    % Set signal to nan for gap and discard events
    signal_interp_curated = process_psa_curate_event("Compute", signal_interp', interp_time, ...
                                                     events, channels, ...
                                                     [], [], {'gap_', 'discard_'})';
    % Fill nan segment with linear interpolation between segment boundaries
    bound_duration = 5;
    interp_time_start = -bound_duration:interp_dt:(interp_time(1)-interp_dt);
    interp_time_end = (interp_time(end)+interp_dt):interp_dt:(interp_time(end)+interp_dt+bound_duration);
    interp_time_ext = [interp_time_start interp_time interp_time_end];

    
    signal_interp_curated_ext = [nan(length(interp_time_start), nb_proc_channels); ...
                                 signal_interp_curated;  ...
                                 nan(length(interp_time_end), nb_proc_channels)];
    signal_interp_filled_ext = fillgaps(signal_interp_curated_ext, 3, 1);

    % Apply spline fitting on downsampled signal filled with linear
    % interpolation in discard and gap segments
    knots_dt = sProcess.options.spline_knot_dt.Value{1};
    knots = interp_time_ext(1):knots_dt:interp_time_ext(end);
    spline_order = 3;
    % spl_mat_sl = spcol(knots, spline_order, interp_time, 'sl');
    spl_mat = spcol(knots, spline_order, interp_time_ext, 'sparse');
    spline_coeffs = spl_mat \ signal_interp_filled_ext;
    signal_spline_fit_ext = spl_mat * spline_coeffs;

    signal_spline_fit = signal_spline_fit_ext(length(interp_time_start)+1:length(signal_interp_curated_ext)-length(interp_time_end), :);

    % Oversample to go back to original time axis
    signal_spline_fit_oversampled = interp1(interp_time', signal_spline_fit, sDataIn.Time', 'linear');

    % figure(); coi=2; hold on; plot(sDataIn.Time, signal(:,coi)', 'k'); plot(interp_time, signal_interp_curated(:,coi)', 'b'); plot(interp_time, signal_spline_fit(:, coi)', 'r'); plot(sDataIn.Time, signal_spline_fit_oversampled(:, coi)', 'g');

    assert(all(all(~isnan(signal_spline_fit_oversampled))));

    if isempty(sProcess.options.spline_channel_names.Value)
        signal_full = sDataIn.F;
        signal_full(chan_mask, :) = signal_spline_fit_oversampled';
    
        % Save time-series data
        sDataOut = db_template('data');
        sDataOut.F            = signal_full;
        sDataOut.Comment      = 'Spline';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = sDataIn.Time;
        sDataOut.History      = sDataIn.History;
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.Events       = events;
        sDataOut.DisplayUnits = sDataIn.DisplayUnits;
        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputs(iInput).iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_spline');
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    else % in-place
        sDataOut = sDataIn;
        chan_names_in = strtrim(str_split(sProcess.options.channel_names.Value, ',;'));
        chan_names_out = strtrim(str_split(sProcess.options.spline_channel_names.Value, ',;'));
        if length(chan_names_in) ~= length(chan_names_out)
            error('Number of destination channels (%d) inconsistent with number of input channels (%d)', ...
                  length(chan_names_out), length(chan_names_in));
        end
        channel_map = containers.Map(chan_names_in, chan_names_out);
        channel_map_idx = containers.Map({chanMat.Channel.Name}, 1:nb_channels);
        processed_chans_idx = find(chan_mask);
        for ichan=1:length(processed_chans_idx)
            chan_idx_in = processed_chans_idx(ichan);
            chan_out_label = channel_map(chanMat.Channel(chan_idx_in).Name);
            chan_out_idx = channel_map_idx(chan_out_label);
            sDataOut.F(chan_out_idx, :) = signal_spline_fit_oversampled(:, ichan);
        end

        for ichan_in=1:length(chan_names_in)
            chan_idx_in = channel_map_idx(chan_names_in{ichan_in});
            chan_out_label = channel_map(chan_names_in{ichan_in});
            chan_idx_out = channel_map_idx(chan_out_label);
            sDataOut.ChannelFlag(chan_idx_out) = sDataOut.ChannelFlag(chan_idx_in);
        end

        OutputFile = file_fullpath(sInputs(iInput).FileName);
        bst_save(OutputFile, sDataOut, 'v7');
    end

    OutputFiles{iInput} = OutputFile;
end
end
