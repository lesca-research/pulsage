function varargout = process_psa_synchronise(varargin )
% PROCESS_EVT_TRANSFER: Synchronise the temporal axis of Input A onto the
% one of Input B.
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
% Authors: Martin Voelker, 2015
%          Raymundo Cassani, 2023
%          Thomas Vincent, 2026

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Synchronise A on B';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Synchronize';
    sProcess.Index       = 681;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    
    %Description of options
    sProcess.options.inputs.Comment = ['Synchronise the time axis of dataset A to the time axis of dataset B.<BR><BR>' ...
        'For synchronization, please choose an<BR>event type which is available in both datasets.<BR><BR>'];
    sProcess.options.inputs.Type    = 'label';
    
    % Source Event name for synchronization (data set A)
    sProcess.options.src.Comment  = 'Sync event name in set A (source): ';
    sProcess.options.src.Type     = 'text';
    sProcess.options.src.Value    = '5';

    % Destination Event name for synchronization (data set B)
    sProcess.options.dest.Comment = 'Sync event name in set B (destination): ';
    sProcess.options.dest.Type    = 'text';
    sProcess.options.dest.Value   = 'E5';

    % Import events from data set B on data set 
    sProcess.options.transfer_events.Comment = 'Event names in set B to import in synchronised set A (comma-separated):';
    sProcess.options.transfer_events.Type    = 'text';
    sProcess.options.transfer_events.Value   = '';  
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB)
    OutputFiles = {};
    
    % Names for sync events
    SyncEventNameA = sProcess.options.src.Value;
    SyncEventNameB = sProcess.options.dest.Value;

    % Load events from datasets A and B
    for iInput = 1 : length(sInputsA)
        sInputA = sInputsA(iInput);
        sInputB = sInputsB(iInput);
        
        % Get file descriptors
        isRawA = strcmpi(sInputA, 'raw');
        if isRawA
            assert(strcmpi(sInputB, 'raw'));
            DataMatA = in_bst_data(sInputA.FileName);
            sFileA = DataMatA.F;
            sEventsA = DataMatA.F.events;
            sFreqA = DataMatA.F.prop.sfreq;
            sFileA.prop.sfreq = sFreqA;
            % Error for epoched CTF and NWB files
            if ismember(sFileA.format, {'CTF', 'NWB'})
                bst_report('Error', sProcess, sInputsA, 'Impossible to process native epoched files. Please import them in database or convert them to continuous.');
            end
            DataMatB = in_bst_data(sInputB.FileName);
            sFileB = DataMatB.F;
            % Error for epoched CTF and NWB files
            if ismember(sFileB.format, {'CTF', 'NWB'})
                bst_report('Error', sProcess, sInputsB, 'Impossible to process native epoched files. Please import them in database or convert them to continuous.');
            end
        else
            DataMatA = in_bst_data(sInputA.FileName);
            sFileA.events = DataMatA.Events;
            sFileA.prop.sfreq = 1 ./ (DataMatA.Time(2) - DataMatA.Time(1));
            sEventsA = DataMatA.Events;
            sFreqA = 1 ./ (DataMatA.Time(2) - DataMatA.Time(1));

            DataMatB = in_bst_data(sInputB.FileName, 'Events', 'Time');
            sFileB.events = DataMatB.Events;
            sFileB.prop.sfreq = 1 ./ (DataMatB.Time(2) - DataMatB.Time(1));
        end

        % Events used for synchronization
        iSyncEventA = find(strcmp({sFileA.events.label}, SyncEventNameA));
        if isempty(iSyncEventA)
            error('Cannot find "%s" in events of Input A', SyncEventNameA);
        end
        iSyncEventB = find(strcmp({sFileB.events.label}, SyncEventNameB));
        if isempty(iSyncEventB)
            error('Cannot find %s in events of Input B', SyncEventNameB);
        end
        nSyncEventA = length(sFileA.events(iSyncEventA).times);
        nSyncEventB = length(sFileB.events(iSyncEventB).times);
    
        % Check for same number of sync events in datasets A and B
        if nSyncEventA ~= nSyncEventB
            error('Cannot process inputs with a different number of synchronization triggers.');
        end
    
        % Calculation of offset between the two datasets
        Offsets  = sFileB.events(iSyncEventB).times - sFileA.events(iSyncEventA).times;
        tOffsetA  = median(Offsets);
        offsetStd = std(Offsets);
        disp(['The population standard deviation of the offsets is ' num2str(offsetStd*1000) ' ms.']);
    
        % Apply offset to dataset A time 
        DataMatA.Time = DataMatA.Time + tOffsetA;
        if isRawA
            DataMatA.F.prop.times = DataMatA.Time;
            if isfield(DataMatA.F, 'epochs') && ~isempty(DataMatA.F.epochs)
                [DataMatA.F.epochs(:).times] = deal(DataMatA.Time);
            end
        end
    
        % Add offset to all events
        for iEvt = 1:length(sEventsA)
            sEventsA(iEvt).times = round((sEventsA(iEvt).times + tOffsetA) .* sFreqA) ./ sFreqA;
        end
        if isRawA
            DataMatA.F.events = sEventsA;
        else
            DataMatA.Events = sEventsA;
        end
    
        %Transfer_events from B to A
        evts_to_transfer = strtrim(str_split(sProcess.options.transfer_events.Value, ',;'));
        if ~isempty(evts_to_transfer)
            found_events = ismember(evts_to_transfer, {sFileB.events.label});
            if ~all(found_events)
                warning('Event(s) to transfer during sync not found: %s', strjoin({evts_to_transfer{~found_events}}, ', '));
            end
            ievts_to_transfer = ismember({sFileB.events.label}, evts_to_transfer);
            sFileA.events = DataMatA.Events;
            sFileA = import_events(sFileA, [], sFileB.events(ievts_to_transfer));
            DataMatA.Events = sFileA.events;
        end
        
        % Add history entry
        DataMatA = bst_history('add', DataMatA, 'timeoffset', sprintf('Added time offset %1.4fs', tOffsetA));
        DataMatA.Comment = [sInputA.Comment ' | sync'];
    
        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputA.iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sync');
        DataMatA.FileName = file_short(OutputFile);
        bst_save(OutputFile, DataMatA, 'v7');
        % Register in database
        db_add_data(sInputA.iStudy, OutputFile, DataMatA);
        OutputFiles{iInput} = OutputFile;
    end
end

