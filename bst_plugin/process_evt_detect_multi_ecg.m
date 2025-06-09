function varargout = process_evt_detect_multi_ecg( varargin )
% PROCESS_EVT_DETECT_MULTI_ECG: Detect heartbeats in a continuous file from several ECG channels,
% and create set of events called "cardiac"

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
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
sProcess.Comment     = 'Detect heartbeats xECG';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Events';
sProcess.Index       = 43;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsDetect#Detection:_Heartbeats';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'raw', 'data'};
sProcess.OutputTypes = {'raw', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Channel name
sProcess.options.channelnames.Comment = 'Channel names (Coma-separated values):';
sProcess.options.channelnames.Type    = 'text';
sProcess.options.channelnames.Value   = '';

% Time window
sProcess.options.timewindow.Comment = 'Time window:';
sProcess.options.timewindow.Type    = 'timewindow';
sProcess.options.timewindow.Value   = [];

% Minimum delay between events
sProcess.options.dt.Comment = 'Torelance interval for peak detection consensus: ';
sProcess.options.dt.Type    = 'value';
sProcess.options.dt.Value   = {0.05, 'ms', 0};

% Event name
sProcess.options.eventname.Comment = 'Event name: ';
sProcess.options.eventname.Type    = 'text';
sProcess.options.eventname.Value   = 'cardiac';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    ChanNames = strtrim(sProcess.options.channelnames.Value);
    if isempty(ChanNames)
        bst_report('Error', sProcess, [], 'Empty channel name.');
        return;
    end
    ChanNames = strtrim(str_split(ChanNames, ',;'));
    if (length(ChanNames) < 1) || isempty(ChanNames{1})
        bst_report('Error', sProcess, [], 'You must enter at least one channel name.');
        return;
    end

    cardiac_event_name = strtrim(sProcess.options.eventname.Value); 

    if length(ChanNames) == 1
        % Process: Detect heartbeats
        OutputFiles = bst_process('CallProcess', 'process_evt_detect_ecg', sInputs, [], ...
                                'channelname', 'ECG I', ...
                                'timewindow',  [], ...
                                'eventname',   cardiac_event_name);
    else
        hb_event_prefix = 'HB_tmp__';
        % Process: Detect heartbeats for every channels
        for ichan=1:length(ChanNames)
            hb_event = [hb_event_prefix num2str(ichan)];
            sFiles = bst_process('CallProcess', 'process_evt_detect_ecg', sInputs, [], ...
                                    'channelname', ChanNames{ichan}, ...
                                    'timewindow',  [], ...
                                    'eventname',   hb_event);
            
        end
        sEvents = in_bst_data(sFiles.FileName, 'Events');
        hb_idx_events = find(~cellfun(@isempty, regexpi({sEvents.Events.label}, [hb_event_prefix '.*'])));
        hb_events = {sEvents.Events(hb_idx_events).label};

        % For peaks that are close together across channels, keep only the
        % first one
        sFiles = bst_process('CallProcess', 'process_evt_multiresp', sFiles, [], ...
                             'responses', strjoin(hb_events, ','), ...
                             'dt',        sProcess.options.dt.Value{1}, ...
                             'action',    1, ...  % Keep only the first event
                             'rename',    0);

        % Merge peaks across channels
        OutputFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
                            'evtnames', strjoin(hb_events, ','), ...
                            'newname',  cardiac_event_name, ...
                            'delete',   1);

    end
end




