function varargout = process_psa_average_ignoring_bad( varargin )

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
% Authors: Aymen Zire (2025-)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Real average pulsatility';
sProcess.FileTag     = 'Average (no bad)';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Average';
sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options

sProcess.options.option_channels.Comment = 'Channels (comma-separated types or names): ';
sProcess.options.option_channels.Type    = 'text';
sProcess.options.option_channels.Value   = '';

sProcess.options.option_win_size.Comment = 'Window size';
sProcess.options.option_win_size.Type    = 'value';
sProcess.options.option_win_size.Value   = {20, 'sec', 2};
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
        % v.1 added
         % masque des bons échantillons 
        sfreq = 1 / diff(sDataIn.Time(1:2));
    nTimes = length(sDataIn.Time);
    good_mask = true(1, nTimes);

    % Liste des labels d'événements à considérer comme "bad"
    bad_event_labels = {'bad', 'bad_segments','BAD','Bad};

    % Appliquer le masque pour chaque événement "bad"
    for iEvt = 1:length(events)
        if any(strcmpi(events(iEvt).label, bad_event_labels))
            for iSeg = 1:size(events(iEvt).times, 2)
                iStart = max(1, round((events(iEvt).times(1, iSeg) - sDataIn.Time(1)) * sfreq) + 1);
                iEnd   = min(nTimes, round((events(iEvt).times(2, iSeg) - sDataIn.Time(1)) * sfreq));
                good_mask(iStart:iEnd) = false;
            end
        end
    end
    
    
    
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iInput).FileName, 'F');
        events = sDataRaw.F.events;
    end
    channels = in_bst_channel(sInputs(iInput).ChannelFile);
    nb_channels = size(channels.Channel, 2);
    if ~isempty(sProcess.options.option_channels.Value)
        idx_chans = channel_find(channels.Channel, sProcess.options.option_channels.Value);
        chan_mask = false(1, nb_channels);
        chan_mask(idx_chans) = 1;
    else
        chan_mask = true(1, nb_channels);
    end

    signal = sDataIn.F(chan_mask, :)';
    win_size = round(sProcess.options.option_win_size.Value{1} / diff(sDataIn.Time(1:2)));
    
    % moyenne glissante ignorant les "bad" OLD signal_mavg = movmean(signal, win_size);
    signal(~good_mask, :) = NaN;
    signal_mavg = movmean(signal, win_size, 'omitnan');
    
    signal_full = sDataIn.F;
    signal_full(chan_mask, :) = signal_mavg';

    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = signal_full;
    sDataOut.Comment      = 'Window-averaged';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time;
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



