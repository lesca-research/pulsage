function varargout = psa_ppl_tcd_autoreg_V1(action, options, arg1, arg2)
%PSA_PPL_TCD_AUTOREG_V1
% Manage an autoregulation pipeline starting from raw TCD + continuous 
% boold pressure data up to autoregulation index computation.
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
%   options = PSA_PPL_TCD_AUTOREG_V1('get_options'); % get default pipeline options
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
%   PSA_PPL_TCD_AUTOREG_V1('save_markings', options, subject_names);
%
%   %% Analysis script
%   options = PSA_PPL_TCD_AUTOREG_V1('get_options');
%   % Define export options again (redundant so this could be factorized)
%   options.export_dir_events = 'path/to/export/events';
%   options.export_dir_channel_flags = 'path/to/export/channel_flags';
%
%   % Customize options:
%   ...
% 
%   % Run the pipeline (and  save user markings):
%   PSA_PPL_TCD_AUTOREG_V1('analyse', options);
%
%   % For a minimal working example see
%   pulsage/script/tcd_autoregulation_pipeline.m [TODO]
%
%% Setup and importation
%
% DEFAULT_OPTIONS = PSA_PPL_TCD_AUTOREG_V1('get_options')
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
        orig_cond = ['origin' get_ppl_tag()];
        files_raw = cellfun(@(s)  nst_get_bst_func_files(s, orig_cond , 'Raw'), ...
                           subject_names, 'UniformOutput', false);
        missing_raw = cellfun(@isempty, files_raw);
        if any(missing_raw)
            warning(sprintf('Missing Raw data files for subjects (will be ignored):\n%s\n', ...
                            strjoin(cellfun(@(s) sprintf(' - %s', s), subject_names(missing_raw), ...
                                    'UniformOutput', false), '\n')));
        end
        files_raw = files_raw(~missing_raw);

        preproc_cond = ['preproc_' get_ppl_tag()];
        files_chan_fix = cellfun(@(s)  nst_get_bst_func_files(s, preproc_cond , 'Raw_Chan_Fix'), ...
                           subject_names, 'UniformOutput', false);
        missing_chan_fix = cellfun(@isempty, files_chan_fix);
        if any(missing_chan_fix)
            warning(sprintf('Missing Raw_Chan_Fix data files for subjects (will be ignored):\n%s\n', ...
                            strjoin(cellfun(@(s) sprintf(' - %s', s), subject_names(missing_chan_fix), ...
                                    'UniformOutput', false), '\n')));
        end
        files_chan_fix = files_chan_fix(~missing_chan_fix);
        
        files_to_sync = [files_raw files_chan_fix];
        
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
pulsatility_indices = struct(); % for table
ar_indices.Subject = '';
for i_cbfi_chan=1:length(CBFi_channel_labels)
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

for iacq=1:length(acq_defs)
    acq_name = acq_defs(iacq).acq_name;
    subject_tag = acq_defs(iacq).subject_tag;

    fprintf('Processing %s\n', acq_name);

    file_raw = nst_get_bst_func_files(acq_name, ['origin' get_ppl_tag()], 'Raw');
    if isempty(file_raw)
        warning('Cannot find "origin/Raw" data for "%s".', acq_name);
        continue;
    elseif iscell(file_raw) && length(file_raw) > 1
        warning('Expect only 1 raw file "origin/Raw", but %d found for "%s".',  ...
            length(file_raw), acq_name);
        continue;
    end

    if iscell(file_raw)
        file_raw = file_raw{1};
    end
        
    data_tmp = in_bst_data(file_raw, 'Events');
    events = data_tmp.Events;

    [~, iStudy] = bst_get('AnyFile', file_raw);
    sChannel = bst_get('ChannelForStudy', iStudy);
    channels = in_bst_channel(sChannel.FileName);
    channel_labels = {channels.Channel.Name};

    cfbi_idx_chans = channel_find(channels.Channel, CBFi_channel_labels);
    if isempty(cfbi_idx_chans)
        error('No CBFi channel found in %s', acq_name);
    end

    if iacq == 1
        i_ar_index = 1;
        ar_indices(i_ar_index).Subject = subject_tag;
    else
        i_ar_index = find(strcmp(subject_tag, {ar_indices.Subject}));  
        if isempty(i_ar_index)
            i_ar_index = length(ar_indices) + 1;
            ar_indices(i_ar_index).Subject = subject_tag;
        end
    end
        
    chan_fix_item = [acq_name '/' preproc_folder '/Raw_Chan_Fix' ];
    if isfield(acq_defs, 'channels_to_swap') && ~isempty(acq_defs(iacq).channels_to_swap)
        [sFile_chan_fix, redone] = nst_run_bst_proc(chan_fix_item, 0, ...
                                                    'process_psa_swap_channels', file_raw, [], ...
                                                    'channelname_1',   acq_defs(iacq).channels_to_swap{1}, ...
                                                    'channelname_2',   acq_defs(iacq).channels_to_swap{2});
    else
        [sFile_chan_fix, redone] = nst_run_bst_proc(chan_fix_item, 0, ...
                                                    'process_duplicate', file_raw, [], ...
                                                    'target', 1, ...  % Duplicate data files
                                                    'tag',    '_copy');
    end

    if options.heart_beats.do
        detect_heart_beats(sFile_chan_fix, ...
                           options.heart_beats.ecg_channel_labels, ...
                           options.heart_beats.event_name);
    end

    if options.do_preproc_only
        continue
    end

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
            raw_cond_item = [acq_name '/' proc_folder '/' condition '_Raw' ];
            [sFile_cond, redone] = nst_run_bst_proc(raw_cond_item, 0, ...
                                                    'process_import_data_event', sFile_chan_fix, [], ...
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

    for ifile=1:length(sFiles_cond)
        % Fill gaps
        sFile_fgaps = sFiles_cond{ifile};
        prefix = prefixes{ifile};
        suffix = suffixes{ifile};

        if options.fill_gaps.do
            for igap=1:length(options.gaps)
                gap_ichannel = find(~cellfun(@isempty, regexp(channel_labels, options.gaps(igap).channel_label)));
                if isempty(gap_ichannel)
                    
                    error('No channel matching %s for gap filling', options.gaps(igap).channel_label);
                end
                fgap_item = [acq_name '/' proc_folder '/' prefix options.gaps(igap).event_label '_filled' ];
                [sFile_fgaps, redone] = nst_run_bst_proc(fgap_item, 0, ...
                                                        'process_psa_fill_gaps', sFile_fgaps, [], ...
                                                        'option_event_name',  options.gaps(igap).event_label, ...
                                                        'option_channels',    channel_labels(gap_ichannel), ...
                                                        'option_ar_win_size', options.fill_gaps.ar_win_size_sec, ...
                                                        'option_ar_order',    options.fill_gaps.ar_order_sec);
            end
        end

        
        mx_item = [acq_name '/' proc_folder '/' prefix 'Mx' ];
        % TODO expose otions for block size and pct filter
        [sFile_mx, redone] = nst_run_bst_proc(mx_item, 0, ...
                                              'process_psa_mx', sFile_fgaps, [], ...
                                              'abp_channel', BP_channel_label, ...
                                              'cbf_channels', channel_labels(cfbi_idx_chans));
        if  options.low_pass_filter.do
            % Low-pass filter
            filter_item = [acq_name '/' proc_folder '/' prefix 'filtered' ];
            [sFile_filtered, redone] = nst_run_bst_proc(filter_item, 0, ...
                'process_bandpass', sFile_fgaps, [], ...
                'sensortypes', 'BP, TCD', ...
                'highpass',    0, ...
                'lowpass',     options.low_pass_filter.low_cutoff, ...
                'tranband',    0, ...
                'attenuation', 'strict', ...  % 60dB
                'ver',         '2019', ...  % 2019
                'mirror',      0, ...
                'overwrite',   0);
        else
            sFile_filtered = sFile_fgaps;
        end

        if options.moving_average.do
            win_avg_comment = [prefix 'mov_average'];
            win_avg_item = [acq_name '/' proc_folder '/' win_avg_comment ];

            [sFile_mov_avg, redone] = nst_run_bst_proc(win_avg_item, 0, ...
                                                        'process_psa_moving_average', sFile_filtered, [], ...
                                                        'option_channels', [BP_channel_label channel_labels(cfbi_idx_chans)], ...
                                                        'option_win_size', options.moving_average_window_sec);
        else
            sFile_mov_avg = sFile_filtered;
        end

        if options.export_preproc.do
            export_fn = fullfile(options.export_dir_preprocessed, ...
                [acq_name '_' win_avg_comment '.tsv']);
            if ~exist(export_fn, 'file')
                export_tmp_fn = fullfile(tempdir, [acq_name '_' win_avg_comment '.tsv']);
                write_log(['Export preprocessed data to tsv for ' acq_name ', ' prefix '\n']);
                export_data(sFile_mov_avg, [], export_tmp_fn, 'ASCII-TSV-HDR-TR');
                copyfile(export_tmp_fn, export_fn);
                delete(export_tmp_fn);
            end
        end
       
        tfa_item = [acq_name '/' proc_folder '/' prefix 'TFA' ];
        [sFile_tfa, redone] = nst_run_bst_proc(tfa_item, 0, ...
                                                'process_psa_tfa', sFile_mov_avg, [], ...
                                                'ref_channel', BP_channel_label, ...
                                                'transfer_channels', channel_labels(cfbi_idx_chans));

        corr_item = [acq_name '/' proc_folder '/' prefix 'corr' ];
        [sFile_corr, redone] = nst_run_bst_proc(corr_item, 0, ...
                                                'process_corr1n', sFile_mov_avg, [], ...
                                                'timewindow',    [], ...
                                                'dest_sensors',  'BP, TCD', ...
                                                'includebad',    0, ...
                                                'timeres',       'none', ...  % None
                                                'avgwinlength',  1, ...
                                                'avgwinoverlap', 50, ...
                                                'scalarprod',    0, ...
                                                'outputmode',    'input');  % separately for each file
        % Load correlation matrix and extract values for CBFi channels
        corr_data = in_bst(sFile_corr);
        corr_mat = process_compress_sym('Expand', corr_data.TF, length(corr_data.RowNames));
        corr_mat = reshape(corr_mat, length(corr_data.RowNames), length(corr_data.RowNames));
        iBP_chan = strcmp(BP_channel_label, corr_data.RowNames);
        MatFlag = in_bst_data(sFile_mov_avg, 'ChannelFlag');

        tfa_data = in_bst(sFile_tfa);
        mx_data = in_bst(sFile_mx);
        for i_tcd_chan=1:length(cfbi_idx_chans)
            tcd_chan_label = channel_labels{cfbi_idx_chans(i_tcd_chan)};
            %corr_tcd_idx_chan = find(~cellfun(@isempty, regexpi(corr_data.RowNames, tcd_chan_label)));
            if iscell(corr_data.RowNames) || isstring(corr_data.RowNames)
                corr_tcd_idx_chan = find(~cellfun(@isempty, regexpi(corr_data.RowNames, tcd_chan_label)));
            else
                warning('corr_data.RowNames is not a cell or string array for subject %s', subject_tag);
                corr_tcd_idx_chan = [];
            end


            if isempty(corr_tcd_idx_chan)
                if MatFlag.ChannelFlag(cfbi_idx_chans(i_tcd_chan)) == 1
                    continue
                    error('TCD channel %s not found in correlation matrix for Mx computation', tcd_chan_label);
                else
                    continue
                end
            end

            mx_index_label = protect_field_label(['AR_Mx_block_' tcd_chan_label]);
            i_mx = strcmp(mx_index_label, mx_data.Description);
            ar_index_label = [mx_index_label suffix];
            if ~isnan(ar_indices(i_ar_index).(ar_index_label))
                error('%s already computed for %s, duplicate condition or duplicate entry in manifest?', ar_index_label, subject_tag);
            end
            % fprintf('Set %s for %s\n', ar_index_label, subject_tag);
            ar_indices(i_ar_index).(ar_index_label) = mx_data.Value(i_mx);

            ar_mx_mova = corr_mat(iBP_chan, corr_tcd_idx_chan);
            ar_index_label = protect_field_label(['AR_Mx_mova_' tcd_chan_label suffix]);
            if ~isnan(ar_indices(i_ar_index).(ar_index_label))
                error('%s already computed for %s, duplicate condition or duplicate entry in manifest?', ar_index_label, subject_tag);
            end
            ar_indices(i_ar_index).(ar_index_label) = ar_mx_mova;

            for itfa_idx=1:length(tfa_index_labels)
                tfa_index_label = protect_field_label([tfa_index_labels{itfa_idx} '_' tcd_chan_label suffix]);
                
                i_tfa = strcmp(tfa_index_label, tfa_data.Description);
                ar_indices(i_ar_index).([tfa_index_label suffix]) = tfa_data.Value(i_tfa);
            end
        end

        if options.pulsatility.do
            pulsatility_item = [acq_name '/' proc_folder '/' prefix 'pulsatility' ];
            
            [sFile_puls, redone] = nst_run_bst_proc(pulsatility_item, 0, ...
                                                    'process_psa_pulsatility', sFile_fgaps, [], ...
                                                    'timewindow',    [], ...
                                                    'channelnames',  strjoin(channel_labels(cfbi_idx_chans), ','), ...
                                                    'heart_beat_event', options.heart_beats.event_name, ...
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
        writetable(pulsatility_table, pulsatility_table_fn, 'WriteRowNames', true, ...
                  'Delimiter', 'tab', ...
                  'FileType', 'text');
        % Génère le graphique à barres groupées à partir du fichier .tsv
        plot_grouped_pulsatility_avg(pulsatility_table_fn);
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

% options.cbf_channel_pattern = 'TCD';

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
ptag = '__psatar_V1';
end

function [files_in, redone_imports] = import_tcd_files(acqs_info, options)
files_in = {};
redone_imports = [];
condition = ['origin' get_ppl_tag()];
for iacq=1:length(acqs_info)
    %% Import data
    tcd_fn = acqs_info(iacq).tcd_fn;
    subject_name = acqs_info(iacq).acq_name;
    discard_start = [];
    if isfield(acqs_info, 'crop_start_time')
        discard_start = acqs_info(iacq).crop_start_time;
    end
    file_raw = nst_get_bst_func_files(subject_name, ...
                                     ['origin' get_ppl_tag()], 'Raw');
    if ~isempty(file_raw)
        files_in = [files_in file_raw];
        redone_imports = [redone_imports 0];
        continue
    else
        redone_imports = [redone_imports 1];
    end
    if endsWith(tcd_fn, 'bst')
        sFiles = bst_process('CallProcess', 'process_import_data_time', [], [], ...
                             'subjectname',  subject_name, ...
                             'condition',    condition, ...
                             'datafile',     {tcd_fn, 'BST-BIN'}, ...
                             'timewindow',   [], ...
                             'split',        0, ...
                             'ignoreshort',  1, ...
                             'channelalign', 1, ...
                             'usectfcomp',   0, ...
                             'usessp',       0, ...
                             'freq',         [], ...
                             'baseline',     []);
    elseif endsWith(tcd_fn, 'tsv')
        % HACK to load TSV file - expect in_fopen_ctf and in_fread_ctf to
        % be overloaded
        sFiles = bst_process('CallProcess', 'process_import_data_time', [], [], ...
            'subjectname',  subject_name, ...
            'condition',    condition, ...
            'datafile',     {tcd_fn, 'CTF'}, ...
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
    sFiles = bst_process('CallProcess', 'process_set_comments', sFiles, [], ...
                         'tag', 'Raw', 'isindex', 0);
    files_in = [files_in sFiles.FileName];
end

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

    for i_cbfi = 1:length(good_cbfi_idx)                   % Boucle sur chaque canal CBFi valide
        chan_idx = good_cbfi_idx(i_cbfi);                  % Index du canal courant
        chan_label = channel_labels{chan_idx};             % Nom du canal

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

        pi_avg = mean(pi_values);                          % Calcule la moyenne

        %fprintf('PI average a.z. : %.4f\n', pi_avg); %TCD left , puis TCD right ET REST first then STAND
        %fprintf('Pulsatilityfile a.z: %s\n', sFile_puls);
        %pi_label = protect_field_label(['Pulsatility_' chan_label suffix]); % Nom du champ de sortie
        %ar_indices(i_ar_index).(pi_label) = pi_avg;        % Stocke la moyenne dans la structure finale
        %pulsatility_indices(i_ar_index).Subject = subject_tag;
        %pulsatility_indices(i_ar_index).(pi_label) = pi_avg;

        % Calcul des statistiques             
        pi_avg = mean(pi_values);

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
        base_label = protect_field_label(['Pulsatility_' chan_label suffix]);                
        % Ajout dans ar_indices
        ar_indices(i_ar_index).([base_label '_avg']) = pi_avg;
        ar_indices(i_ar_index).([base_label '_min']) = pi_min;
        ar_indices(i_ar_index).([base_label '_max']) = pi_max;
        ar_indices(i_ar_index).([base_label '_std']) = pi_std;
        % Ajout dans pulsatility_indices (pour le .tsv)
        pulsatility_indices(i_ar_index).Subject = subject_tag;
        pulsatility_indices(i_ar_index).([base_label '_avg']) = pi_avg;
        pulsatility_indices(i_ar_index).([base_label '_min']) = str_min;
        pulsatility_indices(i_ar_index).([base_label '_max']) = str_max;
        pulsatility_indices(i_ar_index).([base_label '_std']) = pi_std;

        %Two plot
        
       
        base_fig_dir = '/home/lesca-student/ACTIONcR_TCD_autoregulation/data_analysis/pulsatility_figures';

        % Génère un tag unique à partir du nom de fichier complet
        [~, sFileBase, ~] = fileparts(sFile_puls);
        safe_tag = protect_fn_str(sFileBase);  % Nettoie le nom pour l'utiliser dans un dossier
        
        % Utilise un sous-dossier unique pour éviter les collisions
        fig_folder = fullfile(base_fig_dir, [subject_tag '__' safe_tag]);
        
        if ~exist(fig_folder, 'dir')
            mkdir(fig_folder);
        end


        %disp(['Pulsatility file generated: ', sFile_puls]);
        f1 = figure('Visible', 'off');
        boxplot(pi_values);
        title(['Boxplot of PI - ' chan_label suffix]);
        ylabel('Pulsatility Index');
        xlabel('Channel');
        saveas(f1, fullfile(fig_folder, ['boxplot_' chan_label suffix '.png']));
        close(f1);
        
        f2 = figure('Visible', 'off');
        histfit(pi_values, 30);  % Gaussian fit
        title(['Histogram + Fit - ' chan_label suffix]);
        xlabel('Pulsatility Index');
        ylabel('Frequency');
        saveas(f2, fullfile(fig_folder, ['histfit_' chan_label suffix '.png']));
        close(f2);

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
    std_data = T{:, std_labels};

    % https://www.mathworks.com/help/matlab/ref/bar.html
    figure('Visible', 'off');
    bh = bar(avg_data, 'grouped');
    hold on;

    % Calculate number of groups and bars per group
    [nGroups, nBars] = size(avg_data);

    % Find the x positions for each bar
    groupWidth = min(0.8, nBars/(nBars + 1.5));

    errorsbar=true; %afficher ou pas

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
    fig_path = fullfile('/home/lesca-student/ACTIONcR_TCD_autoregulation/data_analysis', 'grouped_pulsatility_avg.png');
    saveas(gcf, fig_path);
    close(gcf);
    disp(['Graphique enregistré ici : ' fig_path]);

end

