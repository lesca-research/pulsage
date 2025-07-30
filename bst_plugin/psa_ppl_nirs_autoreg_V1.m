function varargout = psa_ppl_nirs_autoreg_V1(action, options, arg1, arg2)
%PSA_PPL_NIRS_AUTOREG_V1
% Manage an autoregulation pipeline starting from raw NIRS + continuous 
% blood pressure data up to autoregulation index computation. It was also
% used to get the pulsatility index using ECG and nirs signal.
%
% This pipeline can keep track of user-defined markings outside of brainstorm db 
% such as experiment stages, gap etc. This allows to safely flush all
% brainstorm data while keeping markings.
% 
%% Overview
%
% This function is intended to be called from batch scripts where the user
% can add some custom steps. Here is the workflow:
%
%   %% Importation script
%
%   options = PSA_PPL_NIRS_AUTOREG_V1('get_options'); % get default pipeline options
%
%   % Define options (optional):
%   options.export_dir_events = 'path/to/export/events';
%   options.export_dir_channel_flags = 'path/to/export/channel_flags';
%
%   % Import some data:
%   subject_names = {'subj1', 'subj2'};
%   sFilesRaw = PSA_PPL_TCD_AUTOREG_V1('import', options, {'data1.adicht', 'data2.adicht'}, subject_names);
%   for ifile=1:length(sFilesRaw)
%     % Tweak sFilesRaw{ifile} here, eg import stimulation event.
%   end
%
%   % After the importation script, signals can be manually tagged.
%
%   %% Markings export script
%
%   PSA_PPL_NIRS_AUTOREG_V1('save_markings', options, subject_names);
%
%   %% Analysis script
%   options = PSA_PPL_NIRS_AUTOREG_V1('get_options');
%   % Define export options again (redundant so this could be factorized)
%   options.export_dir_events = 'path/to/export/events';
%   options.export_dir_channel_flags = 'path/to/export/channel_flags';
%
%   % Customize options:
%   ...
% 
%   % Run the pipeline (and  save user markings):
%   PSA_PPL_NIRS_AUTOREG_V1('analyse', options);
%
%   % For a minimal working example see
%   pulsage/script/tcd_autoregulation_pipeline.m [TODO]
%
%% Setup and importation
%
% DEFAULT_OPTIONS = PSA_PPL_NIRS_AUTOREG_V1('get_options')
%     Return default options
%
% FILES_RAW = PSA_PPL_TCD_AUTOREG_V1('import', OPTIONS, ACQS_INFO)
%     Import all nirs files in database and use given subjects (skip if exists).
%
%     Used options:
%        - options.import.redo
%
%     Return:
%         FILES_RAW: brainstorm file pathes to imported data.
%
% FILES_RAW = PSA_PPL_TCD_AUTOREG_V1('import', OPTIONS, LABCHART_FNS, SUBJECT_NAMES)
%     Import all nirs files in database and use given subjects (skip if exists).
%     LABCHART_FNS is a cell array of str.
%     If SUBJECT_NAMES is empty or not given, then use base filename as
%     subject names. If not empty, then it must be a cell array of str with the 
%     same length as LABCHART_FNS.
%
%     Used options:
%        - options.import.redo
%
%     Return:
%         FILES_RAW: brainstorm file pathes to imported data.
%
% FILES_RAW = PSA_PPL_TCD_AUTOREG_V1('save_markings', OPTIONS, SUBJECT_NAMES)
%     Export markings (events and bad channels) for given subjects SUBJECT_NAMES.
%
%     Used options:
%         - options.export_dir_events
%         - options.export_dir_channel_flags
%
%     Return:
%         FILES_RAW: brainstorm file pathes to imported data, for which markings 
%                    were exported.
%
% FILES_RAW = PSA_PPL_TCD_AUTOREG_V1('save_markings', OPTIONS)
%     Export markings (events and bad channels) for all subjects in current protocol.
%
%% Analysis
%
% PSA_PPL_TCD_AUTOREG_V1('analyse', OPTIONS, GROUPS | SUBJECT_NAMES)
%   %
%
global GlobalData;

assert(ischar(action));

switch action
    case 'get_options'
        if nargin > 1
            protocol_name = options;
            assert(ischar(protocol_name));
            varargout{1} = get_options(protocol_name);
        else
            varargout{1} = get_options();
        end
        return;
    case 'import'
        acqs_info = arg1;

        if isempty(options)
            options = get_options();
        end
        create_dirs(options);
        [imported_files, redone] = import_tcd_files(acqs_info, options);
        if options.import.load_markings
            % TODO properly handle timestamp to avoid overwriting when syncing after fresh import 
            %psa_ppl_tcd_autoreg_V1('sync_markings', options, {acqs_info(redone==1).acq_name}, 'load_all');
        end
        varargout{1} = imported_files;
        varargout{2} = redone;
        return;
    case 'import_nirs'
        acqs_info = arg1;
        if isempty(options)
            options = get_options();
        end
        create_dirs(options);
        [imported_files, redone] = import_nirs_files(acqs_info, options);
        %import_markings(imported_files(redone==1), subject_names(redone==1), options);
        % Just in case markings changed in files that were not reimported, save them:
        %export_markings(imported_files(redone==0), subject_names(redone==0), options);
        
        varargout{1} = imported_files;
        varargout{2} = redone;
        return;

    case 'sync_markings'
        if nargin >= 3
            subject_names = arg1;
        else
            % Get all subjects in current protocol (ignore Group_analysis
            % and subject holding full head model)
            if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
                error('Cannot find current protocol');
            end
            sSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
            subject_names = ignore_forged_subjects({sSubjects.Subject.Name});
        end

        if nargin >= 4
            sync_mode = arg2;
        else
            sync_mode = 'sync';
        end

        % TODO adapt for NIRS and Labchart

        %origin saving
        orig_cond = ['origin_labchart' get_ppl_tag()];
        files_raw = cellfun(@(s)  nst_get_bst_func_files(s, orig_cond , 'Raw'), ...
                           subject_names, 'UniformOutput', false);

        missing_raw = cellfun(@isempty, files_raw);
        if any(missing_raw)
            warning(sprintf('Missing Raw data files for subjects (will be ignored):\n%s\n', ...
                            strjoin(cellfun(@(s) sprintf(' - %s', s), subject_names(missing_raw), ...
                                    'UniformOutput', false), '\n')));
        end
        files_raw = files_raw(~missing_raw);

        orig_cond_nirs = ['origin_nirs' get_ppl_tag()];
        files_raw2 = cellfun(@(s)  nst_get_bst_func_files(s, orig_cond_nirs , 'Raw'), ...
                           subject_names, 'UniformOutput', false);

        missing_raw = cellfun(@isempty, files_raw2);
        if any(missing_raw)
            warning(sprintf('Missing Raw data files for subjects (will be ignored):\n%s\n', ...
                            strjoin(cellfun(@(s) sprintf(' - %s', s), subject_names(missing_raw), ...
                                    'UniformOutput', false), '\n')));
        end
        files_raw2 = files_raw2(~missing_raw);


        
        %preproc saving duplicate it for the two preproc folders
        preproc_cond = ['preproc_labchart' get_ppl_tag()];
        files_chan_fix = cellfun(@(s)  nst_get_bst_func_files(s, preproc_cond , 'Events_delete_duplicated_item'), ...
                           subject_names, 'UniformOutput', false);
        missing_chan_fix = cellfun(@isempty, files_chan_fix);
        if any(missing_chan_fix)
            warning(sprintf('Missing Raw_Chan_Fix data files for subjects (will be ignored):\n%s\n', ...
                             strjoin(cellfun(@(s) sprintf(' - %s', s), subject_names(missing_chan_fix), ...
                                    'UniformOutput', false), '\n')));
        end
        files_chan_fix = files_chan_fix(~missing_chan_fix);



        %preproc saving duplicate it for the two preproc folders
        preproc_cond_nirs = ['preproc_nirs' get_ppl_tag()];
        files_chan_fix2 = cellfun(@(s)  nst_get_bst_func_files(s, preproc_cond_nirs , 'Events_uploaded'), ...
                           subject_names, 'UniformOutput', false);
        missing_chan_fix = cellfun(@isempty, files_chan_fix2);
        if any(missing_chan_fix)
            warning(sprintf('Missing Raw_Chan_Fix data files for subjects (will be ignored):\n%s\n', ...
                             strjoin(cellfun(@(s) sprintf(' - %s', s), subject_names(missing_chan_fix), ...
                                    'UniformOutput', false), '\n')));
        end
        files_chan_fix2 = files_chan_fix2(~missing_chan_fix);



        
        files_to_sync = [files_raw files_raw2 files_chan_fix files_chan_fix2];
        
        create_dirs(options);
        sync_markings(files_to_sync, options, sync_mode);
        varargout{1} = files_to_sync;
        return;



    case 'analyse'

    otherwise
        error('Unknown action: %s', action);
end


if isempty(arg1)
    error('Empty input acquisition definitions');
else
    if isstruct(arg1)
        acq_defs = arg1;
        assert(isfield(acq_defs, 'acq_name'));
        assert(isfield(acq_defs, 'subject_tag'));
    else
        acq_defs = repmat(struct('acq_name', '', 'subject_tag', ''), length(arg1), 1);
        for iacq=1:length(arg1)
            acq_defs(iacq).name = arg1{iacq};
            acq_defs(iacq).subject_tag = arg1{iacq};
        end
    end
end

if strcmp(options.save_fig_method, 'export_fig') && ~function_exists('export_fig')
    error('"export_fig" not found. Can be installed from "https://github.com/altmany/export_fig"');
end

create_dirs(options);

protocol_info = bst_get('ProtocolInfo');
if protocol_info.UseDefaultAnat~=1
    error('Protocol should use default anatomy for all subjects');
end

proc_folder = sprintf('proc_%s', get_ppl_tag());
preproc_folder = sprintf('preproc_%s', get_ppl_tag());


BP_channel_label = options.BP_channel_label;
CBFi_channel_labels = options.CBFi_channel_labels;

% Initialise structure to store all AR indices
ar_suffixes = {''};
if ~isempty(options.conditions)     
    for icond=1:length(options.conditions)
        ar_suffixes{icond} = ['_' options.conditions{icond}];
    end
end

tfa_index_labels = process_psa_tfa('tfa_index_labels');
ar_index_labels = [{'AR_Mx_mova', 'AR_Mx_block'} tfa_index_labels];
ar_indices = struct(); 
ar_indices.Subject = '';

for i_cbfi_chan=1:length(CBFi_channel_labels) %itération sur chaque channels NIRS
    cbfi_chan_label = CBFi_channel_labels{i_cbfi_chan};
    for i_ar_index=1:length(ar_index_labels)
        ar_index_label = ar_index_labels{i_ar_index};
        for isuffix=1:length(ar_suffixes)
            ar_chan_index_label = protect_field_label([ar_index_label '_' ...
                                                  cbfi_chan_label ...
                                                  ar_suffixes{isuffix}]);
            ar_indices.(ar_chan_index_label) = nan;
        end
    end
end

pulsatility_indices = struct(); % for table


for iacq=1:length(acq_defs)
    acq_name = acq_defs(iacq).acq_name;
    subject_tag = acq_defs(iacq).subject_tag;

    fprintf('Processing %s\n', acq_name);

    % NIRS Motion correction
    preproc_nirs_folder = ['preproc_nirs' get_ppl_tag()];
    
    % Lacbchart preproc
    origin_labchart_folder = ['origin_labchart' get_ppl_tag()];
    preproc_labchart_folder = ['preproc_labchart' get_ppl_tag()];

    file_raw_labchart = nst_get_bst_func_files(acq_name, origin_labchart_folder, 'Raw');
    if isempty(file_raw_labchart)
        warning('Cannot find "origin_labchart/Raw" data for "%s".', acq_name);
        continue;
    elseif iscell(file_raw_labchart) && length(file_raw_labchart) > 1
        warning('Expect only 1 raw file "origin_labchart/Raw", but %d found for "%s".',  ...
            length(file_raw_labchart), acq_name);
        continue;
    end

    if iscell(file_raw_labchart)
        file_raw_labchart = file_raw_labchart{1};
    end
        
    data_tmp = in_bst_data(file_raw_labchart, 'Events');
    events = data_tmp.Events;

    [~, iStudy] = bst_get('AnyFile', file_raw_labchart);
    sChannel = bst_get('ChannelForStudy', iStudy);
    channels = in_bst_channel(sChannel.FileName);
    channel_labels = {
    'S1D1WL756', 'S3D1WL756', 'S4D1WL756', 'S3D2WL756', 'S4D2WL756', ...
    'S5D2WL756', 'S4D3WL756', 'S5D3WL756', 'S6D3WL756', 'S5D4WL756', ...
    'S6D4WL756', 'S7D4WL756', 'S6D5WL756', 'S7D5WL756', 'S8D5WL756', ...
    'S7D6WL756', 'S8D6WL756', 'S9D6WL756', 'S8D7WL756', 'S9D7WL756', ...
    'S11D7WL756', 'S2D8WL756', 'S10D9WL756', 'S1D1WL853', 'S3D1WL853', ...
    'S4D1WL853', 'S3D2WL853', 'S4D2WL853', 'S5D2WL853', 'S4D3WL853', ...
    'S5D3WL853', 'S6D3WL853', 'S5D4WL853', 'S6D4WL853', 'S7D4WL853', ...
    'S6D5WL853', 'S7D5WL853', 'S8D5WL853', 'S7D6WL853', 'S8D6WL853', ...
    'S9D6WL853', 'S8D7WL853', 'S9D7WL853', 'S11D7WL853', 'S2D8WL853', ...
    'S10D9WL853'
    };
 %not sure if I should run the pipeline for both lamda
% or =opts.psa.CBFi_channel_labels
% 
    


    disp('Tous les canaux disponibles dans le fichier :');
    disp({channels.Channel.Name}');
    cfbi_idx_chans = find(ismember(channel_labels, CBFi_channel_labels));







    chan_fix_item = [acq_name '/' preproc_labchart_folder '/Raw_Fix' ];
    existing_file = nst_get_bst_func_files(acq_name, preproc_labchart_folder, 'Events_delete_duplicated_item');
    
    if ~isempty(existing_file)
        sFile_labchart_hb = existing_file;
        redone = 0;
        fprintf('Skip duplication: Raw_Fix already exists for %s\n', acq_name);
    else
        [sFile_labchart_hb, redone] = nst_run_bst_proc(chan_fix_item, 0, ...
                                                'process_duplicate', file_raw_labchart, [], ...
                                                'target', 1, ...  % Duplicate data files
                                                'tag',    '_copy');
    end
        
        







    if options.heart_beats.do
        detect_heart_beats(sFile_labchart_hb, ...
                           options.heart_beats.ecg_channel_labels, ...
                           options.heart_beats.event_name);
    end

    if options.do_preproc_only
        continue
    end
    
    



    % Merge NIRs + labchart (TODO)
    % === Duplicate to new condition
    process_folder = ['process' get_ppl_tag()];
    target_item = [acq_name '/' process_folder '/Raw'];
    
    
    
    file_raw_nirs = nst_get_bst_func_files(acq_name, ['origin_nirs' get_ppl_tag()], 'Raw');

   
    
    
    
    chan_fix_item = [acq_name '/' preproc_nirs_folder '/Raw_Fix' ];
    existing_file = nst_get_bst_func_files(acq_name, preproc_nirs_folder, 'Events_uploaded');
    
    if ~isempty(existing_file)
        file_nirs_moco = existing_file;
        redone = 0;
        fprintf('Skip duplication: Raw_Fix already exists for %s\n', acq_name);
    else
        [file_nirs_moco, redone] = nst_run_bst_proc(chan_fix_item, 0, ...
                                                'process_duplicate', file_raw_nirs, [], ...
                                                'target', 1, ...  % Duplicate data files
                                                'tag',    '_copy');
    end
        
        
        
        
    
    
    
    
    
    % Renaming events before merging 
    % List of event names to rename as 'Z'
    rename_to_Z_list = {'labchart sync'};  % Add more names here if needed
    % Load existing events
    data_evt_nirs = in_bst_data(file_raw_nirs, 'Events');
    event_labels_nirs = {data_evt_nirs.Events.label};
    
    % Loop over list and rename if present
    for iName = 1:length(rename_to_Z_list)
        src_evt = rename_to_Z_list{iName};
        if any(strcmp(event_labels_nirs, src_evt))
            file_nirs_moco = bst_process('CallProcess', 'process_evt_rename', file_nirs_moco, [], ...
                'src',   src_evt, ...
                'dest',  'Z');
            fprintf('Événement "%s" renommé en "Z" dans %s\n', src_evt, file_nirs_moco.FileName);
            break; 
        end
    end

    % Process: Detect multiple responses and keep one
    evt_delete_double_item = [acq_name '/' 'preproc_labchart__psanar_V1' '/Events_delete_duplicated_item'];
    [file_nirs_proc, ~]= nst_run_bst_proc(evt_delete_double_item,0,'process_evt_multiresp',sFile_labchart_hb, [], ...
                                            'responses', 'CMT oxysoft sync', ...
                                            'dt',        0.001, ...
                                            'action',    1, ...  % Keep only the first event
                                            'rename',    0);
    
    % === Merge events into new copy
    evt_sync_item = [acq_name '/' 'preproc_nirs__psanar_V1' '/Events_uploaded'];
    [file_nirs_proc, ~] = nst_run_bst_proc(evt_sync_item, 0, ...
                                           'process_evt_transfer', ...
                                            sFile_labchart_hb, file_nirs_moco, ...
                                            'src',  'CMT oxysoft sync', ...
                                            'dest', 'Z');

    
    
    continue
    
    sFiles_cond = {};
    prefixes = {};
    suffixes = {};
    if ~isempty(options.conditions)
        for icond=1:length(options.conditions)
            condition = options.conditions{icond};
            idx_event_cond = find(strcmp(condition, {events.label}));
            if isempty(idx_event_cond) || isempty(events(idx_event_cond).times)
                continue
            end
            %condition = ['sync ' options.conditions];
          
            %Importation de rest_RAW et stand_Raw dans proc__psanar_V1
            raw_cond_item = [acq_name '/' proc_folder '/' condition '_Raw' ];
            [sFile_cond, redone] = nst_run_bst_proc(raw_cond_item, 0, ...
                                                    'process_import_data_event', file_nirs_proc, [], ...
                                                    'subjectname',   acq_name, ...
                                                    'condition', '', ...
                                                    'eventname',  condition, ...
                                                    'timewindow',    [], ...
                                                    'epochtime',     [-0.1, 0.3], ...
                                                    'split',         0, ...
                                                    'createcond',    0, ...
                                                    'ignoreshort',   0, ...
                                                    'usectfcomp',    0, ...
                                                    'usessp',        0, ...
                                                    'freq',          [], ...
                                                    'baseline',      [], ...
                                                    'blsensortypes', 'MEG, EEG');
            if isempty(sFile_cond)
                continue
            end
            
            sFiles_cond{end+1} = sFile_cond;
            prefixes{end+1} = [condition '_'];
            suffixes{end+1} = ['_' condition];
        end
    else
        raw_copy_item = [acq_name '/' proc_folder '/' 'Raw' ];
        [sFile_cond, redone] = nst_run_bst_proc(raw_copy_item, 0, ...
                                                'process_duplicate', sFile_chan_fix, [], ...
                                                'target', 1, ...  % Duplicate data files
                                                'tag',    '_copy');

        sFiles_cond = {sFile_cond};
        prefixes = {''};
        suffixes = {''};
    end
    
    i_ar_index = length(pulsatility_indices) + 1;

    for ifile=1:length(sFiles_cond)
        % Fill gaps
        sFile_fgaps = sFiles_cond{ifile};
        prefix = prefixes{ifile};
        suffix = suffixes{ifile};

        
        tag_inversion = [prefix '_INV_nirs'];
        study = bst_get('StudyWithCondition', [acq_name '/' proc_folder]);
        existing_comments = {study.Data.Comment};
        
        if ~any(strcmp(existing_comments, tag_inversion))
            % Process: Opposite values
            sFiles = bst_process('CallProcess', 'process_opposite', sFile_fgaps, [], ...
                                 'sensortypes', 'NIRS', ...
                                 'overwrite', 0);
        
            
            % Process: Set name: _INV_nirs
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                                 'tag', tag_inversion, ...
                                 'isindex', 1);%1 POUR NUMERICAL INDEX
        else
            % Récupérer le fichier existant
            sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
                                 'subjectname', acq_name, ...
                                 'condition',   proc_folder, ...
                                 'tag',         tag_inversion, ...
                                 'includebad',  0, ...
                                 'includeintra', 0);
        end

        tag_motion_correction = [prefix '_motion_correct_nirs'];
        if ~any(strcmp(existing_comments, tag_motion_correction))
            sFiles = bst_process('CallProcess','process_nst_motion_correction', sFiles, [], ...
                                                      'option_event_name', 'movement_artefacts');

            % Process: Set name: _NORM_nirs
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                                 'tag', tag_motion_correction, ...
                                 'isindex', 1);%1 POUR NUMERICAL INDEX
        else
            sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
                                 'subjectname', acq_name, ...
                                 'condition',   proc_folder, ...
                                 'tag',         tag_motion_correction, ...
                                 'includebad',  0, ...
                                 'includeintra', 0);
        end




        tag_normalisation = [prefix '_NORM_nirs'];
        if ~any(strcmp(existing_comments, tag_normalisation))
            % Process: Scale with the mean: [All file]
            sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
                                 'baseline', [], ...
                                 'sensortypes', 'NIRS', ...
                                 'method', 'divmean', ...% Scale with the mean:    x_std = x / &mu;
                                 'overwrite', 0);
            % Process: Set name: _NORM_nirs
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                                 'tag', tag_normalisation, ...
                                 'isindex', 1);%1 POUR NUMERICAL INDEX
        else
            sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
                                 'subjectname', acq_name, ...
                                 'condition',   proc_folder, ...
                                 'tag',         tag_normalisation, ...
                                 'includebad',  0, ...
                                 'includeintra', 0);
        end


        tag_filtre = [prefix '_FILTERED_nirs'];
        if ~any(strcmp(existing_comments, tag_filtre))
            % Process: Band-pass:0.5Hz-5Hz
            sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
                                 'sensortypes', 'NIRS', ...
                                 'highpass',    0.5, ...
                                 'lowpass',     5, ...
                                 'tranband',    0, ...
                                 'attenuation', 'strict', ... % 60dB
                                 'ver',         '2019', ...
                                 'mirror',      0, ...
                                 'overwrite',   0);
            % Process: Set name: _NORM_nirs
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                                 'tag', tag_filtre, ...
                                 'isindex', 1);% Process: Set name: _NORM_nirs
        else
            sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
                                 'subjectname', acq_name, ...
                                 'condition',   proc_folder, ...
                                 'tag',         tag_filtre, ...
                                 'includebad',  0, ...
                                 'includeintra', 0);
        end



        
        if options.pulsatility.do
            pulsatility_item = [acq_name '/' proc_folder '/' prefix 'pulsatility' ];

            [sFile_puls, redone] = nst_run_bst_proc(pulsatility_item, 0, ...
                                                    'process_psa_pulsatility', sFiles, [], ...
                                                    'timewindow',    [], ...
                                                    'channelnames',  strjoin(channel_labels(cfbi_idx_chans), ','), ...
                                                    'heart_beat_event', 'sync cardiac', ...
                                                    'timeres',       'none', ...  % None
                                                    'avgwinlength',  1, ...
                                                    'avgwinoverlap', 50, ...
                                                    'scalarprod',    0, ...
                                                    'outputmode',    'input');  % separately for each file

            %Extrait les statistiques de l'indice de pulsatilité pour chaque canal valide
            % pour chaque canal CBFi valide (non marqué "bad") dans le fichier Brainstorm correspondant.
            % Les résultats sont stockés dans les structures ar_indices et pulsatility_indices pour export ultérieur.
            [ar_indices, pulsatility_indices] = extract_pulsatility_stats_block( sFile_puls, cfbi_idx_chans, channel_labels, ar_indices, ...
                                                    pulsatility_indices, i_ar_index, subject_tag, suffix);
            
        end
    end % end of loop over conditions (eg rest, stand)
end % end of loop over acquisitions

if ~options.do_preproc_only
    ar_table = struct2table(ar_indices, 'AsArray', true); %just added ", 'AsArray', true"
    ar_table_fn = fullfile(options.result_dir, 'autoregulation_indices.tsv');
    write_log(['Save AR index table to ' ar_table_fn '\n']);
    writetable(ar_table, ar_table_fn, 'WriteRowNames', true, ...
              'Delimiter', 'tab', ...
              'FileType', 'text'); 
              
              
    if options.pulsatility.do     
        pulsatility_table = struct2table(pulsatility_indices);
        pulsatility_table_fn = fullfile(options.result_dir, 'pulsatility_indices.tsv');
        write_log(['Save AR index table to ' ar_table_fn '\n']);
        disp(pulsatility_indices);
        writetable(pulsatility_table, pulsatility_table_fn, 'WriteRowNames', true, ...
                  'Delimiter', 'tab', ...
                  'FileType', 'text');
        % Génère le graphique à barres groupées à partir du fichier .tsv
        %plot_grouped_pulsatility_avg(pulsatility_table_fn);
    end
end
end



function detect_heart_beats(sFile, ecg_channel_labels, heart_beat_event_name)
sEvents = in_bst_data(sFile, 'Events');
if ~ismember(heart_beat_event_name, {sEvents.Events.label})
    bst_process('CallProcess', 'process_evt_detect_multi_ecg', sFile, [], ...
                'channelnames', strjoin(ecg_channel_labels, ', '), ...
                'timewindow',   [], ...
                'dt',           0.05, ...
                'eventname',    heart_beat_event_name);
end
end






function label = protect_field_label(label)
label = strrep(label, '.', '_dot_');
end





function [sFiles_flat, tags] = flatten_file_list(sFiles)
tags = {};
sFiles_flat = {};
for isubject=1:length(sFiles)
    files_subject = sFiles{isubject};
    if ~iscell(files_subject)
        files_subject = {files_subject};
    end

    for ifile=1:length(files_subject)
        sFile = in_bst(files_subject{ifile});
        subject_name = fileparts(fileparts(files_subject{ifile}));
        tags = [tags [subject_name '__' sFile.Comment]];
        sFiles_flat = [sFiles_flat files_subject{ifile}];
    end
end

end





function subject_names = ignore_forged_subjects(subject_names)
ignore_subjects = {get_full_head_model_subject_name(), 'Group_analysis'};
subject_names = subject_names(~ismember(subject_names, ignore_subjects));
end




function options = get_options()
% Return default pipeline options

options.redo_all = 0;

options.dry_markings_sync = 0;

options.import.redo = 0;
options.import.ignore_events_in_file = 0;
options.import.load_markings = 1;

options.export_dir_events = ''; % Where to export all events (mirroring), 
                                % espcially those manually defined.
                                % -> will be exported everytime the pipeline
                                %    function is run, or called with the action 
                                %    'export_markings'
                                % -> will be reimported everytime reimportation
                                % options.head_model.surface   is needed.
                              
options.export_dir_channel_flags = ''; % Where to export channel tags (mirroring), 
                                       % espcially those manually tagged.
                                       % -> will be exported everytime the pipeline
                                       %    function is run, or called with the action 
                                       %    'export_markings'
                                       % -> will be reimported everytime reimportation
                                       %    is needed.
options.export_dir_preprocessed = '';

options.do_preproc_only = 0;

options.heart_beats.do = 1;
options.heart_beats.ecg_channel_labels = {'ECG I', 'ECG II', 'ECG III'};
options.heart_beats.event_name = 'cardiac';

options.conditions = {};
options.gaps = struct([]);

options.fill_gaps.do = 1;
options.fill_gaps.ar_win_size_sec = 15;
options.fill_gaps.ar_order_sec = 1;

options.low_pass_filter.do = 1;
options.low_pass_filter.low_cutoff = 0.2; % Hz

options.moving_average.do = 1;
options.moving_average_window_sec = 10;

options.pulsatility.do = 1;

options.deglitch.do = 0;
options.deglitch.redo = 0;
options.deglitch.agrad_std_factor = 2.5;

options.make_figs = 1;
options.save_fig_method = 'saveas'; % 'saveas', 'export_fig'
options.export_fig_dpi = 90;
options.fig_dir = '';
options.fig_background = []; % use default background

end




function ptag = get_ppl_tag()
ptag = '__psanar_V1';
end



function create_dirs(options)
create_dir(options.result_dir);
create_dir(options.fig_dir);
create_dir(options.export_dir_events);
create_dir(options.export_dir_channel_flags);
create_dir(options.export_dir_preprocessed);
end


% TODO: export motion correction tagging to external file
% sRaw = load(file_fullpath(sFile_raw));
% sExport.Events = sRaw.Events(strcmp({sRaw.Events.label}, 'movement_artefacts'));
% export_events(sExport, [], moco_export_fn);

% TODO: export bad channel tagging information

% TODO: plot raw input signals
% fig_bfn = sprintf('%s_%s_signals_raw.png', SubjectName, data_tag);
% fig_fn = protect_fn_str(fullfile(options.fig_dir, fig_bfn ));
% if ~isempty(options.fig_dir) && options.make_figs && options.plot_raw_signals.do && ...
%         (force_redo || options.plot_raw_signals.redo || ~exist(fig_fn, 'file'))
%    plot_signals(sFile_raw, fig_fn, options);
% end


%         if exist(moco_fn, 'file')
            % Load event from pre-saved file
            % TODO: test
%             sFile_in = load(file_fullpath(file_in));
%             [sFile_in, events] = import_events(sFile_in, [], moco_fn, evt_format);   
%         else

function sync_markings(sFiles, options, sync_mode)

io_options.events_to_sync = [options.conditions {options.gaps.event_label, options.heart_beats.event_name}];

if strcmp(sync_mode, 'load_all') 
    bst_process('CallProcess', 'process_evt_delete', sFiles, [], ...
                'eventname', strjoin(io_options.events_to_sync, ','));
end
sync_markings_logic(sFiles, @import_events_wrap, @export_events_wrap, ...
                    io_options, '_events', options.export_dir_events, ...
                    options.dry_markings_sync, sync_mode);

sync_markings_logic(sFiles, @import_channels_wrap, @export_channels_wrap, ...
                    io_options, '_channel_flags', options.export_dir_events, ...
                    options.dry_markings_sync, sync_mode);

end

function export_events_wrap(file_fn, marking_fn, io_options)
nst_bst_export_events(file_fn, marking_fn, io_options.events_to_sync);
end

function import_events_wrap(file_fn, marking_fn, ~)
if exist(marking_fn, 'file') == 2
    nst_bst_import_events(file_fn, marking_fn);
end
end

function import_channels_wrap(file_fn, marking_fn, ~)
if exist(marking_fn, 'file') == 2
    nst_bst_import_channel_flags(file_fn, marking_fn);
end
end

function export_channels_wrap(file_fn, marking_fn, ~)
nst_bst_export_channel_flags(file_fn, marking_fn);
end

function sync_markings_logic(sFiles, import_marking_func, export_marking_func, io_options, marking_fn_suffix, marking_dir, dry, mode)

assert(ismember(mode, {'sync', 'export_all', 'load_all'}));

for iFile=1:length(sFiles)
    file_fn = file_fullpath(sFiles{iFile});
    file_obj = dir(file_fn);
    local_modification_date = file_obj.date;

    [subject_name, condition_label, item_label] = bst_file_info(file_fn);
    marking_fn = fullfile(marking_dir, protect_fn_str([subject_name '--' condition_label '--' item_label ...
                                                       '_' marking_fn_suffix '.mat']));
    if ~exist(marking_fn, 'file')
        export_modification_date = 0;
    else
        sMDate = load(marking_fn, 'modification_date');
        if isempty(fieldnames(sMDate))
            export_modification_date = 0;
        else
            export_modification_date = sMDate.modification_date;
        end
    end
    if strcmp(mode, 'export_all') || (strcmp(mode, 'sync') && local_modification_date > export_modification_date)
        fprintf('Export markings from %s to %s\n', ...
                file_fn, marking_fn);
        if ~dry
            export_marking_func(file_fn, marking_fn, io_options);
            save(marking_fn, 'local_modification_date', '-append');
        end
    elseif strcmp(mode, 'load_all') || (strcmp(mode, 'sync') && local_modification_date < export_modification_date)
        fprintf('Import markings from %s to %s\n', ...
                marking_fn, file_fn);
        if ~dry
            import_marking_func(file_fn, marking_fn, io_options);
        end
    else
        fprintf('No markings sync: Local and export markings have the same date (%s - %s)\n', ...
                marking_fn, file_fn);
    end
end

end

function [subject, condition, comment] = bst_file_info(file_fn)
[aMat, ~] = bst_get('AnyFile', file_fn);
sMat = bst_get('Subject', aMat.BrainStormSubject);
cMat = in_bst_data(file_fn, 'Comment');
subject = sMat.Name;
condition = aMat.Condition{1};
comment = cMat.Comment;
end

function import_markings(sFiles, tags, options)
fprintf('Load events from origin...\n');
io_markings(sFiles, tags, @nst_bst_import_events, @get_events_fn, ...
            'export_dir_events', 1, '', options);
fprintf('Load channel flags from origin...\n');
io_markings(sFiles, tags, @nst_bst_import_channel_flags, @get_channel_flags_fn, ...
            'export_dir_channel_flags', 1, '', options);
end

function export_markings(sFiles, tags, options)
fprintf('Save marking events to origin...\n');
io_markings(sFiles, '_events', ...
            options.export_dir_events, 0, '',  @nst_bst_export_events, [options.conditions {options.gaps.event_label}]);
fprintf('Save channel flags to origin...\n');
io_markings(sFiles, @nst_bst_export_channel_flags, @get_channel_flags_fn, ...
               'export_dir_channel_flags', 0, '', options);
end

function io_markings(sFiles, io_func, get_fn_func, dir_option_name, ...
                      check_exist, msg, options, extra_io_func_arg)
if ~isempty(options.(dir_option_name))
    assert(exist(options.(dir_option_name), 'dir')~=0);
    for iFile=1:length(sFiles)
        sMat = in_bst_data(sFiles{iFile});
        subject_name = sMat.Subject;
        condition_label = sMat.Condition;
        item_label = sMat.Comment;
        marking_fn_prefix(subject_name, condition_label, item_label)
        markings_fn = get_fn_func(marking_fn_prefix,....
                                  options.(dir_option_name));
        if ~check_exist || exist(markings_fn, 'file')==2
            if (nargin >= 9)
                io_func(sFiles{iFile}, markings_fn, extra_io_func_arg);
            else
                io_func(sFiles{iFile}, markings_fn);
            end
            if ~isempty(msg)
                write_log(sprintf([msg '\n'], strrep(markings_fn, '\', '/')));
            end
        end
    end
end
end


% function save_markings_(sFiles, subject_names, export_func, get_fn_func, dir_option_name, msg, options)
% if ~isempty(options.(dir_option_name))
%     assert(exist(options.(dir_option_name), 'dir')~=0);
%     for isubject=1:length(subject_names)
%         markings_fn = get_fn_func(subject_names{isubject},....
%                                   options.(dir_option_name));
%         export_func(sFiles{isubject}, markings_fn);
%         write_log(sprintf('%s to: %s\n', msg, markings_fn));
%     end
% end
% end


function marking_fn = get_marking_fn(bfn, export_dir)
marking_fn = '';
if ~isempty(export_dir)
    assert(exist(export_dir, 'dir')~=0);
    marking_fn = fullfile(export_dir, bfn);
end
end

%% Helper functions

function write_log(msg)
fprintf(msg);
end

function folder = create_dir(folder)
% Create folder if does not exist. 
% Check that folder is not a subfolder of nirstorm sources (encourage good practice
% not to store data in source code folders)

if ~isempty(folder)
    if exist(fullfile(folder, 'nst_install.m'), 'file') || ...
            exist(fullfile(folder, '..', 'nst_install.m'), 'file') || ...
            exist(fullfile(folder, '..', '..', 'nst_install.m'), 'file')
        warning('Processing folder should not be part of nirstorm source folders (%s)', folder);
    end

    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

end

function tag = get_bst_file_tag(fn)

[rr, bfn, ext] = fileparts(fn);
bst_prefixes = {'results_','data_','linkresults_','linkdata_','pdata_','presults_'};
for ipref=1:length(bst_prefixes)
    bfn = replace(bfn, bst_prefixes{ipref}, '');
end

toks = regexp(bfn, '(.*)(?:_\d{6}_\d{4})', 'tokens');
if isempty(toks) || isempty(toks{1})
    tag = bfn;
else
    tag = toks{1}{1};
end

end

function flag = function_exists(func_name)

flag = 1;
try
    eval(func_name);
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        flag = 0;
    end
end

end

function con_str = protect_con_str(con_str)
con_str = strrep(con_str, ' - ', '_minus_');
con_str = strrep(con_str, '-', '_minus_');
con_str = strrep(con_str, ' + ', '_plus_');
con_str = strrep(con_str, '+', '_plus_');
con_str = strrep(con_str, ' * ', '_x_');
con_str = strrep(con_str, '*', '_x_');
con_str = strrep(con_str, ':', '_');
con_str = strrep(con_str, ' ', '_');
end

function sfn = protect_fn_str(s)
sfn = strrep(s, ' | ', '--');
sfn = strrep(s, ' : ', '--');
sfn = strrep(s, ' :', '--');
sfn = strrep(s, ' ', '_');
end


function plot_stat(sFile_ttest, fig_fn, options, show_colbar, save_colbar, sSubjectDefault)
% TODO: set colormap
% TODO: set t-stat cmap boundaries
% TODO: control window size / figure size.
% TODO: contolr displayed scouts
global GlobalData;

hFigSurfData = view_surface_data(sSubjectDefault.Surface(sSubjectDefault.iCortex).FileName, ...
                                 sFile_ttest, 'NIRS', 'NewFigure');
StatThreshOptions = bst_get('StatThreshOptions');
StatThreshOptions.pThreshold = options.GLM_group.contrast_tstat.plot.pvalue_threshold;
StatThreshOptions.Correction = options.GLM_group.contrast_tstat.plot.pvalue_mcc_method;
bst_set('StatThreshOptions', StatThreshOptions);

ColormapInfo = getappdata(hFigSurfData, 'Colormap');
ColormapType = ColormapInfo.Type;
switch ColormapType
    case 'stat2'
        bst_colormaps('SetColormapName', 'stat2', 'cmap_ovun');
    case 'stat1'
        bst_colormaps('SetColormapName', 'stat1', 'cmap_ovun');
end
if ~isempty(options.fig_background)
    bst_figures('SetBackgroundColor', hFigSurfData, options.fig_background);
end
if ~show_colbar
    bst_colormaps('SetDisplayColorbar', 'stat2', 0);
end
view(options.fig_cortex_view);
zoom(hFigSurfData, options.fig_cortex_zoom);

iSurf = getappdata(hFigSurfData, 'iSurface');

panel_surface('SetSurfaceSmooth', hFigSurfData, iSurf, options.fig_cortex_smooth, 0);
% panel_surface('SetSurfaceTransparency', hFigSurfData, iSurf, SurfAlpha);
panel_surface('SetShowSulci', hFigSurfData, iSurf, options.fig_cortex_show_sulci);
panel_scout('SetScoutsOptions', 0, 0, 0, 'none');

%panel_surface('SetDataThreshold',       hFig, iSurf, DataThreshold);
%panel_surface('SetSizeThreshold',       hFig, iSurf, SizeThreshold);

nst_save_figure(fig_fn, options, hFigSurfData);

if save_colbar
    % Save colorbar
    [root, bfn, ext] = fileparts(fig_fn);
    colbar_fig_fn = fullfile(root, [bfn '_colobar' ext]);
    bst_colormaps('SetDisplayColorbar', 'stat2', 1);
    hColorbar = findobj(GlobalData.CurrentFigure.Last, '-depth', 1, 'Tag', 'Colorbar');
    set(hColorbar, 'XColor', [0 0 0]);
    set(hColorbar, 'YColor', [0 0 0]);
    options_colbar = options;
    options_colbar.export_fig_dpi = 500;
    nst_save_figure(colbar_fig_fn, options_colbar, hColorbar);
end
close(hFigSurfData);

end

function [ar_indices, pulsatility_indices] = extract_pulsatility_stats_block( ...
    sFile_puls, cfbi_idx_chans, channel_labels, ar_indices, pulsatility_indices, ...
    i_ar_index, subject_tag, suffix)

    puls_data = in_bst(sFile_puls, [], 1, 1);   % Charge les données de pulsatility en ignore bad segments
    
    %fs = 1 / (puls_data.Time(2) - puls_data.Time(1));  % Calcule la fréquence d'échantillonnage


    
    
    flags_struct = in_bst_data(sFile_puls, 'ChannelFlag'); % Charge les flags des canaux
    
    good_channels = find(flags_struct.ChannelFlag == 1) ;   % prendre indices canaux valides (non marqués "bad") ex : 7,8
    
    good_cbfi_idx = intersect(cfbi_idx_chans, good_channels); % Indices CBFi valides seulement
    %good_cbfi_idx =good_channels;

    for i_cbfi = 1:length(good_cbfi_idx)  % Boucle sur chaque canal CBFi valide
        
        chan_idx = good_cbfi_idx(i_cbfi);                  
        chan_label = channel_labels{chan_idx}; 

        
        % chan_idx = good_cbfi_idx(i_cbfi);
        % % Récupère un label valide pour le canal, peu importe la source
        % if exist('channel_labels', 'var') && numel(channel_labels) >= chan_idx
        %     chan_label = channel_labels{chan_idx};
        % elseif isfield(puls_data, 'RowNames') && numel(puls_data.RowNames) >= chan_idx
        %     chan_label = puls_data.RowNames{chan_idx};
        % else
        %     chan_label = sprintf('Channel_%02d', chan_idx);  % Fallback si aucun label dispo
        % end
        % 
        

        if isfield(puls_data, 'RowNames') && ~isempty(puls_data.RowNames)
            data_idx = find(strcmp(chan_label, puls_data.RowNames)); % Trouve l'index correspondant
        else
            data_idx = chan_idx;                           % Sinon utilise l'index brut
        end

        if isempty(data_idx) || all(isnan(puls_data.F(data_idx, :)))
            warning('No pulsatility data found for channel %s', chan_label); % Si pas de données, sauter
            continue;
        end

        
        signal = puls_data.F(data_idx, :);                 % Signal brut avec NaN
        time_vector = puls_data.Time;                      % Vecteur temps associé
        
        % Exclure les NaN du signal et du temps
        valid_idx = ~isnan(signal);
        pi_values = signal(valid_idx);
        time_vector = time_vector(valid_idx);

        %if isempty(pi_values)
        %    continue;                                      % Ignore si aucune valeur restante
        %end

        pi_avg = median(pi_values);                          % Calcule la moyenne

        %fprintf('PI average a.z. : %.4f\n', pi_avg); %TCD left , puis TCD right ET REST first then STAND
        %fprintf('Pulsatilityfile a.z: %s\n', sFile_puls);
        %pi_label = protect_field_label(['Pulsatility_' chan_label suffix]); % Nom du champ de sortie
        %ar_indices(i_ar_index).(pi_label) = pi_avg;        % Stocke la moyenne dans la structure finale
        %pulsatility_indices(i_ar_index).Subject = subject_tag;
        %pulsatility_indices(i_ar_index).(pi_label) = pi_avg;

        % Calcul des statistiques             
        

        [pi_min, idx_min] = min(pi_values);
        [pi_max, idx_max] = max(pi_values);
        
        % Calcul du temps associé si tu connais ta fréquence d'échantillonnage (ex. 50 Hz)
        %time_vector = (0:length(pi_values)-1) / fs;  % fs = fréquence d'échantillonnage en Hz
        
        time_min = time_vector(idx_min);
        time_max = time_vector(idx_max);
        
        % Formatage texte avec valeur + temps
        str_min = sprintf('%.3f @ %.2fs', pi_min, time_min);
        str_max = sprintf('%.3f @ %.2fs', pi_max, time_max);

        pi_std = std(pi_values);                
        % Construction des noms de champ
        %base_label = protect_field_label(['Pulsatility_' chan_label suffix]); 
        raw_label = ['Pulsatility_' chan_label suffix];
        base_label = matlab.lang.makeValidName(raw_label);  % Rend le nom valide comme champ

        % Ajout dans ar_indices
        ar_indices(i_ar_index).([base_label '_avg']) = pi_avg;
        % ar_indices(i_ar_index).([base_label '_min']) = pi_min;
        % ar_indices(i_ar_index).([base_label '_max']) = pi_max;
        % ar_indices(i_ar_index).([base_label '_std']) = pi_std;
        % Ajout dans pulsatility_indices (pour le .tsv)
        pulsatility_indices(i_ar_index).Subject = subject_tag;
        pulsatility_indices(i_ar_index).([base_label '_avg']) = pi_avg;
        % pulsatility_indices(i_ar_index).([base_label '_min']) = str_min;
        % pulsatility_indices(i_ar_index).([base_label '_max']) = str_max;
        % pulsatility_indices(i_ar_index).([base_label '_std']) = pi_std;

        %Two plot
        
        % changer 
        % base_fig_dir = '/mnt/lesca-data-proc/Project/ACTIONcardioRisk/ACTIONcR_TCD_autoregulation/data_analysis/pulsatility_figures';
        % fig_folder = fullfile(base_fig_dir, subject_tag);
        % if ~exist(fig_folder, 'dir')
        %     mkdir(fig_folder);
        % end
        % 
        % %disp(['Pulsatility file generated: ', sFile_puls]);
        % f1 = figure('Visible', 'off');
        % boxplot(pi_values);
        % title(['Boxplot of PI - ' chan_label suffix]);
        % ylabel('Pulsatility Index');
        % xlabel('Channel');
        % saveas(f1, fullfile(fig_folder, ['boxplot_' chan_label suffix '.png']));
        % close(f1);
        % 
        % f2 = figure('Visible', 'off');
        % histfit(pi_values, 30);  % Gaussian fit
        % title(['Histogram + Fit - ' chan_label suffix]);
        % xlabel('Pulsatility Index');
        % ylabel('Frequency');
        % saveas(f2, fullfile(fig_folder, ['histfit_' chan_label suffix '.png']));
        % close(f2);

    end
end

function plot_grouped_pulsatility_avg(tsv_path)
    [folder_path, file_base, ~] = fileparts(tsv_path);
    T = readtable(tsv_path, 'FileType', 'text', 'Delimiter', '\t');

    % Get averages du tableau et std
    is_avg = endsWith(T.Properties.VariableNames, '_avg');
    avg_data = T{:, is_avg};
    avg_labels = T.Properties.VariableNames(is_avg);
    std_labels = strrep(avg_labels, '_avg', '_std');
    %std_data = T{:, std_labels};

    % https://www.mathworks.com/help/matlab/ref/bar.html
    figure('Visible', 'off');
    bh = bar(avg_data, 'grouped');
    hold on;

    % Calculate number of groups and bars per group
    [nGroups, nBars] = size(avg_data);

    % Find the x positions for each bar
    groupWidth = min(0.8, nBars/(nBars + 1.5));

    errorsbar=false; %afficher ou pas

    if errorsbar
        for i = 1:nBars
            % Calculate center of each bar group
            x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);
            errorbar(x, avg_data(:,i),std_data(:,i), '.', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'CapSize', 3); %errorbar(x, avg_data(:, i), std_data(:, i), 'k.', 'LineWidth', 1);
        end
    end
    title('Pulsatility Index - Moyennes avec écarts-types');
    ylabel('Pulsatility Index moyen');
    xticks(1:nGroups);
    xticklabels(T.Subject);
    xtickangle(45);
    legend(avg_labels, 'Interpreter', 'none', 'Location', 'northeastoutside');
    grid on;
 
    % 
    fig_path = fullfile('/home/lesca-student/ACTIONcR_NIRS_autoregulation/data_analysis', 'grouped_pulsatility_avg.png');
    saveas(gcf, fig_path);
    close(gcf);
    disp(['Graphique enregistré ici : ' fig_path]);

end

%% Pour le import_nirs

function [files_in, redone_imports] = import_nirs_files(acqs_info, options)
files_in = cell(length(acqs_info));
redone_imports = zeros(length(acqs_info));
for iacq=1:length(acqs_info)
    %% Import data
    nirs_fn = acqs_info(iacq).nirs_fn;
    subject_name = acqs_info(iacq).acq_name;
    condition = ['origin_nirs' get_ppl_tag()];
    write_log(['Import nirs data for ' subject_name '\n']);
    [file_in, redone] = nst_run_bst_proc([subject_name '/' condition '/Raw'], options.import.redo, ...
                                           'process_import_data_time', [], [], ...
                                           'subjectname',  subject_name, ...
                                           'condition',    condition, ...
                                           'datafile',     {nirs_fn, 'NIRS-BRS'}, ...
                                           'timewindow',   [], ...
                                           'split',        0, ...
                                           'ignoreshort',  1, ...
                                           'channelalign', 1, ...
                                           'usectfcomp',   0, ...
                                           'usessp',       0, ...
                                           'freq',         [], ...
                                           'baseline',     []);
    % if redone && options.import.ignore_events_in_nirs
    %     data_evt = load(file_fullpath(file_in), 'Events');
    %     sFiles = bst_process('CallProcess', 'process_evt_delete', file_in, [], ...
    %                          'eventname',  strjoin({data_evt.Events.label}, ','));
    % end
    redone_imports(iacq) = redone;
    if redone
        % Create empty event group
        tag_events = db_template('event');
        tag_events(1).label = 'movement_artefacts';
        tag_events(1).times = zeros(2,0);
        tag_events(2).label = 'discard';
        tag_events(2).times = zeros(2,0);
        sFile_in = bst_process('GetInputStruct', file_in);
        process_nst_import_csv_events('import_events', [], sFile_in, tag_events);
    end
    files_in{iacq} = file_in;
end

end


% pour le import
function [files_in, redone_imports] = import_tcd_files(acqs_info, options)
files_in = {};
redone_imports = [];
condition = ['origin_labchart' get_ppl_tag()];
for iacq=1:length(acqs_info)
    %% Import data
    labchart_fn = acqs_info(iacq).labchart_fn;
    subject_name = acqs_info(iacq).acq_name;
    discard_start = [];
    if isfield(acqs_info, 'crop_start_time')
        discard_start = acqs_info(iacq).crop_start_time;
    end
    file_raw = nst_get_bst_func_files(subject_name, ...
                                     ['origin_labchart' get_ppl_tag()], 'Raw');
    if ~isempty(file_raw)
        files_in = [files_in file_raw];
        redone_imports = [redone_imports 0];
        continue
    else
        redone_imports = [redone_imports 1];
    end
    if endsWith(labchart_fn, 'bst')
        sFiles = bst_process('CallProcess', 'process_import_data_time', [], [], ...
                             'subjectname',  subject_name, ...
                             'condition',    condition, ...
                             'datafile',     {labchart_fn, 'BST-BIN'}, ...
                             'timewindow',   [], ...
                             'split',        0, ...
                             'ignoreshort',  1, ...
                             'channelalign', 1, ...
                             'usectfcomp',   0, ...
                             'usessp',       0, ...
                             'freq',         [], ...
                             'baseline',     []);
    elseif endsWith(labchart_fn, 'tsv')
        % HACK to load TSV file - expect in_fopen_ctf and in_fread_ctf to
        % be overloaded
        sFiles = bst_process('CallProcess', 'process_import_data_time', [], [], ...
            'subjectname',  subject_name, ...
            'condition',    condition, ...
            'datafile',     {labchart_fn, 'CTF'}, ...
            'timewindow',   [], ...
            'split',        0, ...
            'ignoreshort',  1, ...
            'channelalign', 1, ...
            'usectfcomp',   0, ...
            'usessp',       0, ...
            'freq',         [], ...
            'baseline',     []);
    else
        error('Unsupported file format');
    end
    warning off MATLAB:load:variableNotFound
    data_tmp = in_bst_data(sFiles.FileName, 'Events');
    warning on MATLAB:load:variableNotFound
    if isempty(data_tmp.Events) || ~any(ismember({options.gaps.event_label}, {data_tmp.Events.label}))
        tag_events = db_template('event');
        ievt = 1;
        for igap=1:length(options.gaps)
            tag_events(ievt).label = options.gaps(igap).event_label;
            tag_events(ievt).times = zeros(2,0);
            ievt = ievt + 1;
        end
        for icondition=1:length(options.conditions)
            tag_events(ievt).label = options.conditions{icondition};
            tag_events(ievt).times = zeros(2,0);
            ievt = ievt + 1;
        end
        % sFile_in = bst_process('GetInputStruct', sFiles);
        process_nst_import_csv_events('import_events', [], sFiles, tag_events);
    end

    if ~isempty(discard_start) && ~isnan(discard_start)
        DataMatOrig = in_bst_data(sFiles.FileName);
        fprintf('Crop %s from %d sec\n',  subject_name, discard_start);
        time_window = [discard_start DataMatOrig.Time(end)];
        sFiles = bst_process('CallProcess', 'process_extract_time', sFiles, [], ...
                            'timewindow', time_window, ...
                            'overwrite',  1);
        % DataMat = in_bst_data(sFiles.FileName);
        % DataMat.AllUnits = DataMatOrig.AllUnits;
        % DataMat.Timestamp = DataMatOrig.Timestamp;
        % bst_save(file_fullpath(sFiles.FileName), DataMat, 'v7');
    end
    sFiles = bst_process('CallProcess', 'process_set_comments', sFiles, [], ...
                         'tag', 'Raw', 'isindex', 0);
    files_in = [files_in sFiles.FileName];
end

end
