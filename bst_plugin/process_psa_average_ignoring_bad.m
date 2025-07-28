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
sProcess.Index       = 1304; %0: not shown, >0: defines place in the list of processes
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
     
    
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iInput).FileName, 'F');
        events = sDataRaw.F.events;
    end
 % v.1 added
     % masque des bons échantillons 
    sfreq = 1 / diff(sDataIn.Time(1:2));
    nTimes = length(sDataIn.Time);
    good_mask = true(1, nTimes);

     % Liste des labels d'événements à considérer comme "bad"
    bad_event_labels = {'bad', 'bad_segments','BAD','Bad'};

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



    channels = in_bst_channel(sInputs(iInput).ChannelFile);
    nb_channels = size(channels.Channel, 2);
    if ~isempty(sProcess.options.option_channels.Value)
        idx_chans = channel_find(channels.Channel, sProcess.options.option_channels.Value);
        chan_mask = false(1, nb_channels);
        chan_mask(idx_chans) = 1;
    else
        chan_mask = true(1, nb_channels);
    end

    signal = sDataIn.F(chan_mask, :)'; % signal : time x channel


        % === Moyenne AVEC les bad ===
    signal_orig = signal;  % version complète du signal
    mean_with_bad = mean(signal_orig, 1);  % moyenne classique
    
    % === Moyenne SANS les bad ===
    signal_bad_masked = signal;  % copie à modifier
    signal_bad_masked(~good_mask, :) = NaN;  % NaN sur les segments "bad"
    signal_avg = nanmean(signal_bad_masked, 1);  % moyenne sans bad, 1 x channel
    
    % === Format texte pour commentaire du fichier ===
    mean_with_bad_str = num2str(mean_with_bad(1), '%.2f');
    mean_no_bad_str   = num2str(signal_avg(1), '%.2f');

    fprintf('[%s] Mean with bad: %.4f | Mean without bad: %.4f\n', sInputs(iInput).FileName, mean_with_bad(1), signal_avg(1));

    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = signal_avg';
    sDataOut.Comment = ['Mean (bad incl: ' mean_with_bad_str ', excl: ' mean_no_bad_str ')'];
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = [1];
    sDataOut.History      = sDataIn.History;
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
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



