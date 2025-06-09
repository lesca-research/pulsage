function varargout = process_psa_fill_gaps( varargin )

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
sProcess.Comment     = 'Fill Gaps';
sProcess.FileTag     = 'fill gap';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Artifacts';
sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options
sProcess.options.option_event_name.Comment = 'Name of event group defining gaps: ';
sProcess.options.option_event_name.Type    = 'text';
sProcess.options.option_event_name.Value   = 'Gap';

sProcess.options.option_channels.Comment = 'Channels (comma-separated types or names): ';
sProcess.options.option_channels.Type    = 'text';
sProcess.options.option_channels.Value   = '';

sProcess.options.option_ar_win_size.Comment = 'Window size for AR estimation';
sProcess.options.option_ar_win_size.Type    = 'value';
sProcess.options.option_ar_win_size.Value   = {15, 'sec', 2};   

sProcess.options.option_ar_order.Comment = 'AR order (~duration of one period)';
sProcess.options.option_ar_order.Type    = 'value';
sProcess.options.option_ar_order.Value   = {1, 'sec', 2};
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};

if isempty(which('fillgaps'))
    bst_error('Curve Fitting Toolbox not available');
    return 
end
% elseif isempty(which('csaps'))
%     bst_error(['Curve Fitting Toolbox OK but function csaps not found.<BR>' ...
%                'Try refreshing matlab cache using command: rehash toolboxcache']);
%     return
% end

% Get selected events
event_name =  strtrim(sProcess.options.option_event_name.Value);

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
    
    event = [];
    ievt_fgap = [];
    for ievt=1:length(events)
        if strcmp(events(ievt).label, event_name)
            event = events(ievt);
            ievt_fgap = ievt;
            break;
        end
    end
    if isempty(event)
        warning(['Event "' event_name '" does not exist in file.']);
    end

    if isempty(event) || isempty(event.times) % no marked event
        signal_filled = sDataIn.F';
    else
        dt = diff(sDataIn.Time(1:2));
        ar_win_len = round(sProcess.options.option_ar_win_size.Value{1} / dt);
        ar_order = round(sProcess.options.option_ar_order.Value{1} / dt);
        channels = in_bst_channel(sInputs(iInput).ChannelFile);
        nb_channels = size(channels.Channel, 2);
        
        if ~isempty(sProcess.options.option_channels.Value)
            idx_chans = channel_find(channels.Channel, sProcess.options.option_channels.Value);
            main_chan_mask = false(1, nb_channels);
            main_chan_mask(idx_chans) = 1;
        else
            main_chan_mask = true(1, nb_channels);
        end

        signal = sDataIn.F';
        samples = time_to_sample_idx(event.times, sDataIn.Time);
        to_fill = false(size(signal));
        for ievt=1:size(event.times, 2)
            if isempty(event.channels) || isempty(event.channels{ievt}) % apply to all channels
                chan_mask = main_chan_mask;
            else
                % ignore main channel filter and use channel of event
                chan_mask = strcmp(event.channels{ievt}, {channels.Channel.Name});
            end
            to_fill(samples(1, ievt):samples(2, ievt), chan_mask) = 1;
        end
        signal_filled = Compute(signal, ar_win_len, to_fill, ar_order);
    end
    
    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = signal_filled';
    sDataOut.Comment      = 'Gap-filled';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    if ~isempty(ievt_fgap)
        sDataOut.Events       = events([1:(ievt_fgap-1) (ievt_fgap+1):length(events)]);
    else
        sDataOut.Events       = events;
    end
    sDataOut.DisplayUnits = sDataIn.DisplayUnits;
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_motion_corr');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    OutputFiles{iInput} = OutputFile;
end
end


%% ===== Compute =====
function [sig_filled] = Compute(sig, AR_win_len, fill_mask, AR_order) %#ok<DEFNU>
%
% Args
%    - sig: matrix of double, size: time x nb_channels
%        signals to be corrected
%    - AR_win_len: integer (nb of samples)
%        Window size of the left and right signal chunks around the gap
%        on which the AR model is fitted.
%    - AR_order: integer (nb of samples)
%        Order of the AR model, analog to its memory.
%        For periodic signals, it is advice to choose the length of one
%        period.
%    - fill_mask: matrix of boolean, size: time x nb_channels
%        boolean mask indicating gaps to be filled
%
% Output:
%    - sig_filled: matrix of double, size: time x nb_channels
%        Filled signals
%
extra_nans = isnan(sig);
extra_nans(fill_mask) = 0;
if any(extra_nans)
    warning('Nans in input signal will be temporary filled with zeros before filling the given gaps');
    sig(extra_nans) = 0;
end
sig(fill_mask) = nan;
sig_filled = fillgaps(sig, AR_win_len, AR_order);
sig_filled(extra_nans) = nan;
end

function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end
