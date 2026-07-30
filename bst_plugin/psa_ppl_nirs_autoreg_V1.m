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
    case 'create_template'
        create_template(arg1);
        return;    
    case 'import'
        acqs_info = arg1;

        if isempty(options)
            options = get_options();
        end
        create_dirs(options);
        [imported_files, redone] = import_data_files(acqs_info, options);
        varargout{1} = imported_files;
        varargout{2} = redone;
        return;
    case 'save_markings'
        acqs_info = arg1;
        ppl_tag = get_ppl_tag();
        files_to_sync = {};
        for iacq=1:length(acqs_info)
            subject_name = acqs_info(iacq).subject_name;
            visit = acqs_info(iacq).visit_label;
            condition_nirs = [visit '_origin_nirs' ppl_tag];
            file_nirs_raw = nst_get_bst_func_files(subject_name, condition_nirs, 'Raw');
            if isempty(file_nirs_raw)
                warning('Missing NIRS Raw data files for %s - %s\n', ...
                        subject_name, visit);
            else
                files_to_sync{end+1} = file_nirs_raw{1};
            end
            condition_aux = [visit '_origin_aux' ppl_tag];
            file_aux_raw = nst_get_bst_func_files(subject_name, condition_aux, 'Raw');
            if isempty(file_aux_raw)
                warning('Missing AUX Raw data files for %s - %s\n', ...
                        subject_name, visit);
            else
                files_to_sync{end+1} = file_aux_raw{1};
            end
        end
        create_dirs(options);
        sync_markings(files_to_sync, options, 'export_all');
        varargout{1} = files_to_sync;
        return;

     case 'load_markings'
        acqs_info = arg1;
        ppl_tag = get_ppl_tag();
        files_to_sync = {};
        for iacq=1:length(acqs_info)
            subject_name = acqs_info(iacq).subject_name;
            visit = acqs_info(iacq).visit_label;
            condition_nirs = [visit '_origin_nirs' ppl_tag];
            file_nirs_raw = nst_get_bst_func_files(subject_name, condition_nirs, 'Raw');
            if isempty(file_nirs_raw)
                error('Missing NIRS Raw data files for %s - %s\n', ...
                        subject_name, visit);
            else
                files_to_sync{end+1} = file_nirs_raw{1};
            end
            condition_aux = [visit '_origin_aux' ppl_tag];
            file_aux_raw = nst_get_bst_func_files(subject_name, condition_aux, 'Raw');
            if isempty(file_aux_raw)
                error('Missing AUX Raw data files for %s - %s\n', ...
                        subject_name, visit);
            else
                files_to_sync{end+1} = file_aux_raw{1};
            end
        end
        create_dirs(options);
        sync_markings(files_to_sync, options, 'load_all');
        varargout{1} = files_to_sync;
        return;
    case 'detect_heart_beats'
        detect_heart_beats(arg1, options)
        return;
    case 'analyse'

    otherwise
        error('Unknown action: %s', action);
end

if isempty(arg1)
    error('Empty input acquisition definitions');
else
    acq_defs = arg1;
    assert(isfield(acq_defs, 'visit'));
    assert(isfield(acq_defs, 'subject_name'));
end

if strcmp(options.save_fig_method, 'export_fig') && ~function_exists('export_fig')
    error('"export_fig" not found. Can be installed from "https://github.com/altmany/export_fig"');
end

create_dirs(options);

protocol_info = bst_get('ProtocolInfo');
if protocol_info.UseDefaultAnat~=1
    error('Protocol should use default anatomy for all subjects');
end

if isempty(options.BP_channel_label)
    error('options.BP_channel_label is undefined');
end
BP_channel_label = options.BP_channel_label;

all_features = dictionary();
ppl_tag = get_ppl_tag();
for iacq=1:length(acq_defs)
    subject_name = acq_defs(iacq).subject_name;
    visit = acq_defs(iacq).visit;
    
    fprintf('Processing %s - %s\n', subject_name, visit);

    orig_cond_nirs = [visit '_origin_nirs' ppl_tag];
    file_nirs_raw = nst_get_bst_func_files(subject_name, orig_cond_nirs , 'Raw');
    file_nirs_raw = file_nirs_raw{1};
    orig_cond_aux = [visit '_origin_aux' ppl_tag];
    file_aux_raw = nst_get_bst_func_files(subject_name, orig_cond_aux , 'Raw');
    file_aux_raw = file_aux_raw{1};

    % TODO QC mark cardiac events that are too close or too spread apart

    % Remove duplicate sync triggers
    nirs_time = in_bst(file_nirs_raw, 'Time');
    bst_process('CallProcess', 'process_evt_multiresp', ...
                 file_nirs_raw, [], ...
                'responses', options.sync_event.nirs, ...
                'dt',        diff(nirs_time(1:2)), ...
                'action',    1, ...  % Keep only the first event
                'rename',    0);
    nirs_mat = in_bst_data(file_nirs_raw, 'Events');
    nirs_events = nirs_mat.Events;

    aux_time = in_bst(file_aux_raw, 'Time');
    bst_process('CallProcess', 'process_evt_multiresp', ...
                 file_aux_raw, [], ...
                'responses', options.sync_event.aux, ...
                'dt',        diff(aux_time(1:2)), ...
                'action',    1, ...  % Keep only the first event
                'rename',    0);
    aux_mat = in_bst_data(file_aux_raw, 'Events');
    aux_events = aux_mat.Events;

    % Insure that conditions have been tagged
    if ~isempty(options.conditions)
        if isempty([aux_events(ismember({aux_events.label}, options.conditions)).times])
            warning('Condition (%s) not tagged in %s/%s/Raw', strjoin(options.conditions, ', '), ...
                subject_name, orig_cond_aux);
        elseif ~isempty([nirs_events(ismember({nirs_events.label}, options.conditions)).times])
            error('Condition should not be tagged in both nirs and aux');
        end
    end
    % Synchronise NIRS onto AUX
    if ismember(options.sync_event.nirs, {nirs_events.label}) && ismember(options.sync_event.aux, {aux_events.label})
        file_nirs_sync = nst_run_bst_proc(fmt_item(visit, 'procs', 'nirs', 'Raw_sync_aux'), 0, ...
                                         'process_psa_synchronise',  file_nirs_raw, file_aux_raw, ...
                                         'src',  options.sync_event.nirs, ...
                                         'dest', options.sync_event.aux, ...
                                         'transfer_events', strjoin([options.sync_event.aux options.conditions options.heart_beats.event_name 'discard_BP' 'discard_cardiac' 'extrasystole'], ','));
        nirs_mat_tmp = in_bst_data(file_nirs_sync, 'Time');
        aux_mat_tmp = in_bst_data(file_aux_raw, 'Time');
        if 0
            fprintf('Time bounds of nirs after sync: %1.2f - %1.2f\n', ...
                    nirs_mat_tmp.Time(1), nirs_mat_tmp.Time(end));
            fprintf('Time bounds of aux: %1.2f - %1.2f\n', ...
                    aux_mat_tmp.Time(1), aux_mat_tmp.Time(end));
        end
        common_tstart = max(nirs_mat_tmp.Time(1), aux_mat_tmp.Time(1));
        common_tend = min(nirs_mat_tmp.Time(end), aux_mat_tmp.Time(end));
        % fprintf('Common time bounds: %1.2f - %1.2f\n', common_tstart, common_tend);
        file_nirs_aligned = nst_run_bst_proc(fmt_item(visit, 'procs', 'nirs', 'Raw_aligned_aux'), 0, ...
                                             'process_import_data_time', file_nirs_sync, [], ...
                                             'subjectname',  subject_name, ...
                                             'condition',     '', ...
                                             'timewindow',    [common_tstart, common_tend], ...
                                             'split',         0, ...
                                             'ignoreshort',   0, ...
                                             'usectfcomp',    0, ...
                                             'usessp',        0, ...
                                             'freq',          [], ...
                                             'baseline',      [], ...
                                             'blsensortypes', 'MEG, EEG');
         file_aux_aligned = nst_run_bst_proc(fmt_item(visit, 'procs', 'aux', 'Raw_aligned_nirs'), 0, ...
                                             'process_import_data_time', file_aux_raw, [], ...
                                             'subjectname',  subject_name, ...
                                             'condition',     '', ...
                                             'timewindow',    [common_tstart, common_tend], ...
                                             'split',         0, ...
                                             'ignoreshort',   0, ...
                                             'usectfcomp',    0, ...
                                             'usessp',        0, ...
                                             'freq',          [], ...
                                             'baseline',      [], ...
                                             'blsensortypes', 'MEG, EEG');
        sync_ok = true;
    else
        warning('No sync event available. Autoregulation will not be computed.');
        sync_ok = false;
        file_nirs_aligned = file_nirs_raw;
    end
    data_nirs_tmp = in_bst_data(file_nirs_aligned, 'Events', 'Time');
    nirs_events = data_nirs_tmp.Events;

    % NIRS: Crop trailing segments tagged with event "discard"
    idiscard_evt = strcmp('discard', {nirs_events.label});
    file_nirs_discard = file_nirs_aligned;
    if any(idiscard_evt) && ~isempty(nirs_events(idiscard_evt).times)
        discard_evt = nirs_events(idiscard_evt);
        data_time = in_bst_data(file_nirs_aligned, 'Time');
        time_window = [0 data_time.Time(end)];
        if discard_evt.times(1,1)==data_time.Time(1)
            time_window(1) = discard_evt.times(2,1);
            discard_evt.times(:,1) = [];
        end
        if ~isempty(discard_evt.times) && discard_evt.times(2,end) == data_time.Time(end)
            time_window(2) = discard_evt.times(1,end);
            discard_evt.times(:,end) = [];
        end
        if time_window(1) ~= 0 || time_window(2) ~= data_time.Time(end)
            file_nirs_discard = nst_run_bst_proc(fmt_item(visit, 'procs', 'nirs', 'Discard'), 0, ...
                                                  'process_extract_time', file_nirs_aligned, [], ...
                                                  'timewindow', time_window, ...
                                                  'overwrite',  0);
            bst_process('CallProcess', 'process_evt_delete', file_nirs_discard, [], ...
                        'eventname', 'discard');
            if ~isempty(discard_evt)
                events_tmp = db_template('event');
                events_tmp(idiscard_evt) = discard_evt;
                sInput_tmp.FileType = 'data';
                sInput_tmp.FileName = sFile_cropped_discard;
                process_nst_import_csv_events('import_events', [], sInput_tmp, events_tmp);
            end
        end
    end
    % NIRS Deglitching
    file_nirs_deglitched = nst_run_bst_proc(fmt_item(visit, 'procs', 'nirs', 'Deglitched'), 0, ...
                                            'process_nst_deglitch', file_nirs_discard, [], ...
                                            'factor_std_grad', options.deglitch.agrad_std_factor);
    
    % NIRS Motion corrections
    file_nirs_moco = nst_run_bst_proc(fmt_item(visit, 'procs', 'nirs', 'Moco'), 0, ...
                                     'process_nst_motion_correction', file_nirs_deglitched, [], ...
                                     'option_event_name', 'movement_artefacts');
    % fprintf('Split NIRS by condition\n');
    [sFiles_nirs_conds, conds] = split_by_condition(file_nirs_moco, options.conditions, 'Moco', options);

    if sync_ok
        % BP gap filling
        [~, iStudy] = bst_get('AnyFile', file_aux_raw);
        sChannel = bst_get('ChannelForStudy', iStudy);
        channels = in_bst_channel(sChannel.FileName);
        channel_labels = {channels.Channel.Name};
        gap_ichannel = find(~cellfun(@isempty, regexp(channel_labels, BP_channel_label)));
        if isempty(gap_ichannel)
            error('No channel matching %s for gap filling', BP_channel_label);
        end
        fgap_item = fmt_item(visit, 'procs', 'aux', 'BP_gap_filled');
        [file_aux_fgaps, ~] = nst_run_bst_proc(fgap_item, 0, ...
                                            'process_psa_fill_gaps', file_aux_aligned, [], ...
                                            'option_event_name',  'BP_gap', ...
                                            'option_channels',    BP_channel_label, ...
                                            'option_ar_win_size', options.fill_gaps.ar_win_size_sec, ...
                                            'option_ar_order',    options.fill_gaps.ar_order_sec);            
        % fprintf('Split AUX by condition\n');
        [sFiles_aux_conds, aux_conds] = split_by_condition(file_aux_fgaps, options.conditions, 'BP_gap_filled', options);
    
        if length(conds) ~= length(aux_conds)
            error('Inconsitent conditions between nirs (%s) and aux (%s)', ...
                  none_str(strjoin(conds, ', ')), none_str(strjoin(aux_conds, ', ')));
        end
        % Check consitency of cardiac events in NIRS and AUX
        for icond=1:length(sFiles_nirs_conds)
            data_nirs_tmp = in_bst_data(sFiles_nirs_conds{icond}, 'Events', 'Time');
            nirs_events = data_nirs_tmp.Events;
            
            data_aux_tmp = in_bst_data(sFiles_aux_conds{icond}, 'Events', 'Time');
            aux_events = data_aux_tmp.Events;
    
            dt_nirs = diff(data_nirs_tmp.Time(1:2));
            if any(abs(data_nirs_tmp.Time([1 end]) - data_aux_tmp.Time([1 end])) > dt_nirs) 
                warning('Different %s time bounds btwn NIRS (%1.2f-%1.2f) and AUX (%1.2f-%1.2f)', ...
                        conds{icond}, data_nirs_tmp.Time(1), data_nirs_tmp.Time(end), ...
                        data_aux_tmp.Time(1), data_aux_tmp.Time(end));
            end
            icardiac_nirs = find(strcmp({nirs_events.label}, options.heart_beats.event_name));
            icardiac_aux = find(strcmp({aux_events.label}, options.heart_beats.event_name));
            nirs_cardiac = nirs_events(icardiac_nirs).times();
            aux_cardiac = aux_events(icardiac_aux).times();
            if any(abs(nirs_cardiac([1 end]) - aux_cardiac([1 end])) > dt_nirs)
                 warning('Different %s cardiac event bounds btwn NIRS (%1.2f-%1.2f) and AUX (%1.2f-%1.2f)', ...
                         conds{icond}, nirs_cardiac(1), nirs_cardiac(end), ...
                          aux_cardiac(1), aux_cardiac(end));
            end
        end

    end

    % Pulsatility NIRS
    for icond=1:length(sFiles_nirs_conds)
        condition = conds{icond};
        sFile_nirs = sFiles_nirs_conds{icond};
        [pfile, pfeature_means, pfeature_labels] = pulsatility_pipeline(sFile_nirs, condition, options.heart_beats.event_name);
        for ifeature=1:length(pfeature_labels)
            key = {{subject_name, visit, [pfeature_labels{ifeature} '_' condition]}};
            all_features(key) = pfeature_means(ifeature);
        end
    end

    if ~sync_ok
        continue % don't do autoregulation
    end

    for icond=1:length(sFiles_nirs_conds)
        condition = conds{icond};
        sFile_nirs = sFiles_nirs_conds{icond};
        sFile_aux = sFiles_aux_conds{icond};
        
        data_aux_tmp = in_bst_data(sFile_aux, 'Events', 'Time', 'ChannelFlag');
        aux_events = data_aux_tmp.Events;
        idiscard_evt = find(strcmp('discard_BP', {aux_events.label}));
        if ~isempty(idiscard_evt) && sum(diff(aux_events(idiscard_evt).times)) / diff(data_aux_tmp.Time([1 end])) > 0.5
            warning('Ignore AUX BP %s signal because more than 50% is discarded', ...
                    condition);
            continue;
        end
        ichan_BP = find(~cellfun(@isempty, regexp(channel_labels, BP_channel_label)));
        if isempty(ichan_BP) || data_aux_tmp.ChannelFlag(ichan_BP) ~= 1
            warning('No BP available for %s - %s', subject_name, condition)
            continue;
        end
        % TODO expose option for pct filter
        mx_item = fmt_item(visit, 'procs', 'nirs', [condition '_Mx']);
        [sFile_mx, redone] = nst_run_bst_proc(mx_item, 0, ...
                                              'process_psa_mx2', sFile_nirs, sFile_aux, ...
                                              'abp_channel', BP_channel_label);
        [~, iStudy] = bst_get('AnyFile', sFile_mx);
        sChannel = bst_get('ChannelForStudy', iStudy);
        channels = in_bst_channel(sChannel.FileName);
        channels_selection = strjoin({channels.Channel(endsWith({channels.Channel.Name}, 'WL853')).Name}, ',');
        mx_mean = channel_average(sFile_mx, 0.02, channels_selection, 5);
        key = {{subject_name, visit, ['Mx_' condition]}};
        all_features(key) = mx_mean;
    end
end % end of loop over acquisitions

%% Export to Excel
feature_table = dict_to_table(all_features, {'subject_id', 'visit'});
outputFile = fullfile(options.result_dir, 'hemodynamics.xlsx');
writetable(feature_table, outputFile, 'Sheet', 'Hemodynamics_Features');
fprintf('Wrote %d rows and %d columns to %s\n', ...
        height(feature_table), width(feature_table), outputFile);
disp(feature_table);
end % end main pipeline function


function table_out = dict_to_table(dict, var_names)

allKeys     = keys(dict);     % cell array of 1×4 cells
allMeasures = values(dict);   % numeric column vector

% Unpack each heterogeneous key into its components
n_entries = numel(allKeys);

% Pre-allocate index columns as cells (since types may vary)
nb_vars = length(var_names);
columns = cell(nb_vars, n_entries);
labelCol = strings(n_entries, 1);   % measure_label is always a string

for i_entry = 1:n_entries
    thisKey = allKeys{i_entry};
    if length(thisKey)-1 ~= nb_vars
        error('Entry size in dict inconstistent with number of given variables');
    end
    for  icol=1:nb_vars
        columns{icol, i_entry} = thisKey{icol};
    end
    labelCol(i_entry) = string(thisKey{icol+1});
end

% Convert each index column to its most natural MATLAB type
columns_normed = cell(n_entries, nb_vars+2);
for icol=1:nb_vars
    col = normalizeColumn(columns(icol, :));
    for i_entry=1:n_entries
        columns_normed{i_entry, icol} = col{i_entry};
    end
end
for i_entry=1:n_entries
    columns_normed{i_entry, icol+1} = labelCol{i_entry};
    columns_normed{i_entry, icol+2} = allMeasures(i_entry);
end
% Long-format table
table_long = cell2table(columns_normed, 'VariableNames', [var_names {'measure_label','measure_value'}]);

% Pivot to wide format: one column per measure label
table_out = unstack(table_long, 'measure_value', 'measure_label');
end

%% Helper: collapse a cell column to a homogeneous type when possible
function out = normalizeColumn(c)
    % If every element is numeric scalar, return a numeric column.
    % If every element is char/string, return a string column.
    % Otherwise, leave it as a cell column (mixed types).
    
    isAllNumeric = all(cellfun(@(x) isnumeric(x) && isscalar(x), c));
    isAllText    = all(cellfun(@(x) ischar(x) || isstring(x), c));
    
    if isAllNumeric
        out = cell2mat(c);
    elseif isAllText
        out = string(c);
    else
        out = c;   % keep heterogeneous as cell
    end
end

function s = none_str(s)
if isempty(s)
    s = 'none';
end
end

function [subject_name, proc_folder] = parse_bst_file(bst_file)
tokens = strsplit(bst_file, '/');
proc_folder = tokens{2};
subject_name = tokens{1};
end

function [sFile_features_avg, feature_means, feature_labels] = pulsatility_pipeline(file_moco, tag, cardiac_event_label)
[subject_name, proc_folder] = parse_bst_file(file_moco);
if isempty(tag)
    prefix = '';
else
    prefix = [tag '_'];
end
item_prefix = [subject_name '/' proc_folder '/' prefix];
[sFile_baseline, redone] = nst_run_bst_proc([item_prefix 'baseline_LP'], 0, ...
                                             'process_bandpass', file_moco, [], ...
                                             'sensortypes', 'NIRS', ...
                                             'highpass',    0, ...
                                             'lowpass',     0.6, ...
                                             'tranband',    0, ...
                                             'attenuation', 'strict', ... % 60dB
                                             'ver',         '2019', ...
                                             'mirror',      0, ...
                                             'overwrite',   0);
if isempty(tag)
    condition_tag = 'lp_normed';
else
    condition_tag = tag;
end
division_code = ['Data = (DataA ./ DataB - 1)*100;' 10 'Condition = ''' condition_tag ''';'];
[sFile_normed, ~] = nst_run_bst_proc([item_prefix 'LP_normalised'], 0, ...
                                          'process_matlab_eval2', file_moco, sFile_baseline, ...
                                          'matlab', division_code, ...
                                          'sensortypes', 'NIRS');
channels = in_bst_channel(bst_get('ChannelFileForStudy', sFile_normed));
separations = process_nst_separations('Compute', channels.Channel);

feature_labels = {'pulse_amplitude', 'tt_pCBF', 'pCBF', 'heart_rate'};
feature_means = zeros(length(feature_labels), 1);
sFile_features_avg = {};
for ifeature=1:length(feature_labels)
    feature = feature_labels{ifeature};
    [sFile_feat, ~] = nst_run_bst_proc([item_prefix feature], 0, ...
                                        'process_psa_pulsatility', sFile_normed, [], ...
                                        'timewindow',    [], ...
                                        'channelnames',  'NIRS', ...
                                        'heart_beat_event', cardiac_event_label, ...
                                        'use_discard_channel_events',   1, ...
                                        'discard_global_events', 'discard_cardiac, discard, extrasystole', ...
                                        'flip_sign', 1, ...
                                        'function', feature);
    [sFile_feat_avg, redone] = nst_run_bst_proc([item_prefix feature '_avg'], 0, ...
                                                'process_psa_average_time_robust', sFile_feat, [], ...
                                                'pct_filter_range', [25, 75], ...
                                                'avg_func',  'mean', ...  % Arithmetic average:  mean(x)
                                                'overwrite',        0);
    sFile_features_avg{ifeature} = sFile_feat_avg;
    channels_selection = strjoin({channels.Channel(endsWith({channels.Channel.Name}, 'WL853')).Name}, ',');
    % % Compute average across all channels
    % feature_avg_data = in_bst(sFile_feat_avg);
    % idx_chans = channel_find(channels.Channel, 'NIRS');
    % mask_good_chans = false(size(feature_avg_data.ChannelFlag, 1), 1);
    % mask_good_chans(idx_chans) = 1;
    % % TODO expose wavelength of interest
    % % TODO expose separation threshold
    % mask_good_chans(separations < 0.02) = 0;
    % mask_good_chans(~endsWith({channels.Channel.Name}, 'WL853')) = 0;
    % mask_good_chans(feature_avg_data.ChannelFlag~=1) = 0;
    % %TODO filter by SCI
    % feature_data = feature_avg_data.F(mask_good_chans', 1);
    % nb_channels = size(feature_data, 1);
    % low = max(1, round(0.05 * nb_channels));
    % high = round(.95 * nb_channels);
    % [sorted_data, sort_idx] = sort(feature_data);
    % feature_mean = mean(sorted_data(low:high));
    feature_means(ifeature) = channel_average(sFile_feat_avg, 0.02, channels_selection, 5);
end
end

function feature_mean = channel_average(data_file, separation_max, channels_selection, pct_filter)

% Compute average across all channels
channels = in_bst_channel(bst_get('ChannelFileForStudy', data_file));
separations = process_nst_separations('Compute', channels.Channel);
data = in_bst(data_file);
idx_chans = channel_find(channels.Channel, channels_selection);
mask_good_chans = false(size(data.ChannelFlag, 1), 1);
mask_good_chans(idx_chans) = 1;
% TODO expose wavelength of interest
% TODO expose separation threshold
mask_good_chans(separations < separation_max) = 0;
% mask_good_chans(~endsWith({channels.Channel.Name}, 'WL853')) = 0;
mask_good_chans(data.ChannelFlag~=1) = 0;
%TODO filter by SCI
feature_data = data.F(mask_good_chans', 1);
nb_channels = size(feature_data, 1);
pct_filter = pct_filter/100;
low = max(1, round(pct_filter * nb_channels));
high = round((1-pct_filter) * nb_channels);
[sorted_data, sort_idx] = sort(feature_data);
feature_mean = mean(sorted_data(low:high));
end

function [output_files, output_conditions] = split_by_condition(data_file, conditions, tag, options)
data_mat = in_bst_data(data_file, 'Events');
events = data_mat.Events;
[subject_name, proc_folder] = parse_bst_file(data_file);
if ~isempty(conditions)
    output_files = {};
    output_conditions = {};
    for icond=1:length(conditions)
        condition = conditions{icond};
        idx_event_cond = find(strcmp(condition, {events.label}));
        if isempty(idx_event_cond) || isempty(events(idx_event_cond).times)
            continue
        end
        if 0
            fprintf('Split by %s: t=[%1.2f-%1.2f], dur=%1.2f\n', condition, ...
                    events(idx_event_cond).times(1), ...
                    events(idx_event_cond).times(2), ...
                    events(idx_event_cond).times(2)-events(idx_event_cond).times(1));
        end
        idx_event_cardiac = find(strcmp({events.label}, ...
                                 options.heart_beats.event_name));
        if ~isempty(idx_event_cardiac)
            cardiac_times = events(idx_event_cardiac).times;
            nb_cardiac = sum((cardiac_times >= events(idx_event_cond).times(1)) & ...
                             (cardiac_times <= events(idx_event_cond).times(2)));
            if 0
                fprintf('Nb of cardiac events in %s before split: %d\n', condition, ...
                        nb_cardiac);
            end
        end
        raw_cond_item = [subject_name '/' proc_folder '/' condition '_' tag];
        [sFile_cond, ~] = nst_run_bst_proc(raw_cond_item, 0, ...
                                                'process_import_data_event', data_file, [], ...
                                                'subjectname',   subject_name, ...
                                                'condition', proc_folder, ...
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
       
        data_mat_tmp = in_bst_data(sFile_cond, 'Events', 'Time');
        events_after = data_mat_tmp.Events;
        if 0
        fprintf('Time bound of extracted %s signal: t=[%1.2f-%1.2f], dur=%1.2f\n',...
                condition, data_mat_tmp.Time(1), data_mat_tmp.Time(end), ...
                data_mat_tmp.Time(end) - data_mat_tmp.Time(1));
        end
        idx_event_cardiac = find(strcmp({events_after.label}, ...
                                 options.heart_beats.event_name));
        if 0 && ~isempty(idx_event_cardiac)
            fprintf('Nb of cardiac events in %s after split: %d\n', condition, ...
                    length(events_after(idx_event_cardiac).times));
        end

        output_files{end+1} = sFile_cond;
        output_conditions{end+1} = condition;
    end
else
    [data_file_dup, ~] = nst_run_bst_proc([subjet_name '/' proc_folder '/Raw' ], 0, ...
                                            'process_duplicate', data_file, [], ...
                                            'target', 1, ...  % Duplicate data files
                                            'tag', '_copy');

    output_files = {data_file_dup};
    output_conditions = {''};
end

end

function item=fmt_item(visit, proc, signal_type, tag)
    item = sprintf('%s_%s_%s%s/%s', visit, proc, signal_type,...
                   get_ppl_tag(), tag);
end

function create_template(nirs_fn)
fm_subject = get_full_head_model_subject_name();
file_raw_fm = nst_get_bst_func_files(fm_subject, ['origin' get_ppl_tag()], 'Raw');
if isempty(file_raw_fm)
    if endsWith(nirs_fn, '.snirf')
        format = 'NIRS-SNIRF';
    else
        format = 'NIRS-BRS';
    end
    sFile_in = bst_process('CallProcess', 'process_import_data_time', [], [], ...
                            'subjectname',  fm_subject, ...
                            'condition',    ['origin' get_ppl_tag()], ...
                            'datafile',     {nirs_fn, format}, ...
                            'timewindow',   [], ...
                            'split',        0, ...
                            'ignoreshort',  1, ...
                            'channelalign', 1, ...
                            'usectfcomp',   0, ...
                            'usessp',       0, ...
                            'freq',         [], ...
                            'baseline',     []);
    sFile_in = bst_process('CallProcess', 'process_set_comment', sFile_in, [], ...
                            'tag',     'Raw', ...
                            'isindex', 0);
    
    optodes_coords_fn = fullfile(fileparts(nirs_fn), 'template__nspst_V1_optodes.tsv');
    if exist(optodes_coords_fn, "file") == 2
        fprintf('Set optode coordinates from file: %s\n', optodes_coords_fn);
        bst_process('CallProcess', 'process_channel_addloc', sFile_in, [], ...
                    'channelfile', {optodes_coords_fn, 'BIDS-SCANRAS-MM'}, ...
                    'usedefault',  '', ...  %
                    'fixunits',    1, ...
                    'vox2ras',     0, ...
                    'mrifile',     {'', ''}, ...
                    'fiducials',   []);
        manual_fix = 0;
    else
        fprintf('No optode coordinates file found (template__nspst_V1_optodes.tsv)');
        manual_fix = 1;
    end
    fprintf('Template created.\n');
    if manual_fix
        fprintf('You can manually adjust the optode montage:\n'); 
        fprintf('Right-click on template channels > MRI Registration > Edit...\n');
        fprintf('Once the montage is adjusted, it is advised to export the optode coordinates.\n')
    end
    fprintf('Next step: Import all data\n');
else
    fprintf('Template already created. Remove it before calling this function to create it again.');
end

end



function add_default_nirs_events(sFile, options)
if ~isstruct(sFile)
    sFile = bst_process('GetInputStruct', sFile);
end
warning off MATLAB:load:variableNotFound
data_tmp = in_bst_data(sFile.FileName, 'Events');
warning on MATLAB:load:variableNotFound
if ~isempty(data_tmp.Events)
    data_event_labels = {data_tmp.Events.label};
else
    data_event_labels = {};
end
tag_events = db_template('event');
 
ievt = 1;
if ~ismember('movement_artefacts', data_event_labels)
    tag_events(ievt).label = 'movement_artefacts';
    tag_events(ievt).times = zeros(2,0);
    ievt = ievt + 1;
end
if ~ismember('discard', data_event_labels)
    tag_events(ievt).label = 'discard';
    tag_events(ievt).times = zeros(2,0);
    ievt = ievt + 1;
end
if ievt > 1
    process_nst_import_csv_events('import_events', [], sFile, ...
                                  tag_events);
end
end

function add_default_aux_events(sFile, options)
if ~isstruct(sFile)
    sFile = bst_process('GetInputStruct', sFile);
end
warning off MATLAB:load:variableNotFound
data_tmp = in_bst_data(sFile.FileName, 'Events');
warning on MATLAB:load:variableNotFound
if ~isempty(data_tmp.Events)
    data_event_labels = {data_tmp.Events.label};
else
    data_event_labels = {};
end
tag_events = db_template('event');

ievt = 1;
% Segments where BP is artefactual and can be reconstructed < 2 sec.
if ~ismember('BP_gap', data_event_labels)
    tag_events(ievt).label = 'BP_gap';
    tag_events(ievt).times = zeros(2,0);
    ievt = ievt + 1;
end
% Segments where BP is artefactual and are too long (> 2 sec) to be reconstructed
if ~ismember('discard_BP', data_event_labels)
    tag_events(1).label = 'discard_BP';
    tag_events(1).times = zeros(2,0);
    ievt = ievt + 1;
end
% Segments where heart beats cannot be tagged
if ~ismember('discard_cardiac', data_event_labels)
    tag_events(1).label = 'discard_cardiac';
    tag_events(1).times = zeros(2,0);
    ievt = ievt + 1;
end

if ~ismember('extrasystole', data_event_labels)
    tag_events(ievt).label = 'extrasystole';
    tag_events(ievt).times = zeros(2,0);
    ievt = ievt + 1;
end
for icondition=1:length(options.conditions)
    if ~ismember(options.conditions{icondition}, data_event_labels)
        tag_events(ievt).label = options.conditions{icondition};
        tag_events(ievt).times = zeros(2,0);
        ievt = ievt + 1;
    end
end
% sFile_in = bst_process('GetInputStruct', sFiles);
if ievt > 1
    process_nst_import_csv_events('import_events', [], sFile,...
                                  tag_events);
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
options.heart_beats.good_ecg_event_name = 'good_ecg';

options.BP_channel_label = '';

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

io_options.events_to_sync = [options.conditions {'cardiac', 'movement_artefacts', 'discard' , 'BP_gap', 'discard_BP', 'discard_cardiac', 'extrasystole'}];

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
if strcmp(sync_mode, 'load_all') 
    sInputs = bst_process('GetInputStruct', sFiles);
    for ifile=1:length(sFiles)
        if ~isempty(strfind(sInputs(ifile).FileName, 'nirs'))
            add_default_nirs_events(sInputs(ifile), options);
        elseif ~isempty(strfind(sInputs(ifile).FileName, 'aux'))
            add_default_aux_events(sInputs(ifile), options);
        else
            error('Cannot recognise signal type %s', sInputs(ifile).FileName);
        end
    end
end
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

assert(ismember(mode, {'export_all', 'load_all'}));

for iFile=1:length(sFiles)
    file_fn = file_fullpath(sFiles{iFile});
    [subject_name, condition_label, item_label] = bst_file_info(file_fn);
    marking_fn = fullfile(marking_dir, protect_fn_str([subject_name '--' condition_label '--' item_label ...
                                                       '_' marking_fn_suffix '.mat']));
    if strcmp(mode, 'export_all') 
        fprintf('Export markings from %s to %s\n', ...
                file_fn, marking_fn);
        if ~dry
            export_marking_func(file_fn, marking_fn, io_options);
        end
    elseif strcmp(mode, 'load_all') 
        fprintf('Import markings from %s to %s\n', marking_fn, file_fn);
        if ~dry
            import_marking_func(file_fn, marking_fn, io_options);
        end
    else
        error('Bad mode %s\n', mode);
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
            options.export_dir_events, 0, '',  @nst_bst_export_events, [options.conditions {options.gaps.event_label, options.discard_segments.event_label}]);
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

function fm_subject = get_full_head_model_subject_name()
    fm_subject = ['template' get_ppl_tag()];
end

% pour le import
function [files_in, redone_imports] = import_data_files(acqs_info, options)

% Export template channel coordinates. 
% They will be applied to all imported files
fm_subject = get_full_head_model_subject_name();
sStudy = bst_get('StudyWithCondition', [fm_subject '/origin' get_ppl_tag()]);
template_ChannelFile = file_fullpath(sStudy.Channel.FileName);
tmp_template_channel_fn = [tempname '.tsv'];
export_channel(template_ChannelFile, tmp_template_channel_fn, 'BIDS-NIRS-ALS-MM', 0);

files_in = {};
redone_imports = [];

ppl_tag = get_ppl_tag();
for iacq=1:length(acqs_info)
    %% Import data
    subject_name = acqs_info(iacq).subject_name;
    visit = acqs_info(iacq).visit_label;
    aux_fn = acqs_info(iacq).labchart_fn;
    nirs_fn = acqs_info(iacq).nirs_fn;
    
    condition_nirs = [visit '_origin_nirs' ppl_tag];
    file_nirs_raw = nst_get_bst_func_files(subject_name, condition_nirs, 'Raw');
    if isempty(file_nirs_raw)
        [file_nirs_raw, redone] = nst_run_bst_proc([subject_name '/' condition_nirs '/Raw'], options.import.redo, ...
                                       'process_import_data_time', [], [], ...
                                       'subjectname',  subject_name, ...
                                       'condition',    condition_nirs, ...
                                       'datafile',     {nirs_fn, 'NIRS-BRS'}, ...
                                       'timewindow',   [], ...
                                       'split',        0, ...
                                       'ignoreshort',  1, ...
                                       'channelalign', 1, ...
                                       'usectfcomp',   0, ...
                                       'usessp',       0, ...
                                       'freq',         [], ...
                                       'baseline',     []);
        % Insure events loaded from NIRS have prefix CMT
        rename_matching_event(file_nirs_raw, @(s) (~startsWith(s, 'CMT')), ...
                              @(s) (['CMT ' s]));
        % warning off MATLAB:load:variableNotFound
        % data_tmp = in_bst_data(file_nirs_raw, 'Events');
        % warning on MATLAB:load:variableNotFound
        % to_rename = {};
        % renamed = {};
        % for ievt=1:length(data_tmp.Events)
        %     if ~startsWith(data_tmp.Events(ievt).label, 'CMT')
        %         to_rename{end+1} = data_tmp.Events(ievt).label;
        %         renamed{end+1} = ['CMT ' data_tmp.Events(ievt).label];
        %     end
        % end
        % if ~isempty(to_rename)
        %     bst_process('CallProcess', 'process_evt_rename', file_nirs_raw, [], ...
        %         'src',   strjoin(to_rename, ', '), 'dest',  strjoin(renamed, ','));
        % end
        merge_duplicate_events(file_nirs_raw);
        % TODO apply proper channel coordinates?
    else
        file_nirs_raw = file_nirs_raw{1};
    end

    condition_aux = [visit '_origin_aux' ppl_tag];
    file_aux_raw = nst_get_bst_func_files(subject_name, condition_aux, 'Raw');
    if isempty(file_aux_raw)
        [file_aux_raw, redone] = nst_run_bst_proc([subject_name '/' condition_aux '/Raw'], options.import.redo, ...
                                                'process_import_data_time', [], [], ...
                                                'subjectname',  subject_name, ...
                                                'condition',    condition_aux, ...
                                                'datafile',     {aux_fn, 'BST-BIN'}, ...
                                                'timewindow',   [], ...
                                                'split',        0, ...
                                                'ignoreshort',  1, ...
                                                'channelalign', 1, ...
                                                'usectfcomp',   0, ...
                                                'usessp',       0, ...
                                                'freq',         [], ...
                                                'baseline',     []);
    else
        file_aux_raw = file_aux_raw{1};
    end
    merge_duplicate_events(file_aux_raw);

    add_default_nirs_events(file_nirs_raw, options)
    add_default_aux_events(file_aux_raw, options);
end

detect_heart_beats(acqs_info, options);
end

function merge_duplicate_events(file_raw)
% Merge events with duplicate group names
warning off MATLAB:load:variableNotFound
data_tmp = in_bst_data(file_raw, 'Events');
warning on MATLAB:load:variableNotFound
if (~isfield(data_tmp, 'Events') || ...
    isempty(data_tmp.Events) || ...
    length(unique({data_tmp.Events.label})) == length(data_tmp.Events))
    return;
end
merged_evt = dictionary();
for ievt=1:length(data_tmp.Events)
    evt = data_tmp.Events(ievt).label;
    if ~isConfigured(merged_evt) || ~isKey(merged_evt, evt)
        merged_evt(evt) = {data_tmp.Events(ievt).times};
    else
        previous_events = merged_evt(evt);
        merged_evt(evt) = {[previous_events{1} data_tmp.Events(ievt).times]};
    end
end

new_events = db_template('event');
evts = keys(merged_evt);
for ievt=1:length(evts)
    evt = evts{ievt};
    evt_times = merged_evt(evt);
    new_events(ievt).label = evt;
    new_events(ievt).times = evt_times{1};
    new_events(ievt).epochs = ones(1, size(evt_times{1}, 2));
end

bst_process('CallProcess', 'process_evt_delete', file_raw, [], ...
            'eventname', strjoin({data_tmp.Events.label}, ','));
process_nst_import_csv_events('import_events', [],...
                               bst_process('GetInputStruct', file_raw), new_events);

end

function rename_matching_event(data_file, match_func, rename_func)
data_tmp = in_bst_data(data_file, 'Events');
to_rename = {};
renamed = {};
for ievt=1:length(data_tmp.Events)
    if match_func(data_tmp.Events(ievt).label)
        to_rename{end+1} = data_tmp.Events(ievt).label;
        renamed{end+1} = rename_func(data_tmp.Events(ievt).label);
    end
end
if ~isempty(to_rename)
    bst_process('CallProcess', 'process_evt_rename', data_file, [], ...
                'src',   strjoin(to_rename, ', '), 'dest',  strjoin(renamed, ','));
end

end

function detect_heart_beats(acqs_info, options)
ppl_tag = get_ppl_tag();
for iacq=1:length(acqs_info)
    %% Import data
    subject_name = acqs_info(iacq).subject_name;
    visit = acqs_info(iacq).visit_label;
    
    condition_nirs = [visit '_origin_nirs' ppl_tag];
    file_nirs_raw = nst_get_bst_func_files(subject_name, condition_nirs, 'Raw');
    
    if ~isempty(file_nirs_raw)
        file_nirs_raw = file_nirs_raw{1};
        data_nirs_tmp = in_bst_data(file_nirs_raw, 'Events', 'Time');
        nirs_events = data_nirs_tmp.Events;
        if (isfield(acqs_info(iacq), 'cardiac_channel_ref') && ...
            ~(isempty(acqs_info(iacq).cardiac_channel_ref)))
            detect_heart_beats_nirs(file_nirs_raw, acqs_info(iacq).cardiac_channel_ref);
        elseif ~ismember(options.sync_event.nirs, {nirs_events.label})
            error('No sync event in NIRS and cardiac reference channel not defined for %s - %s.', ...
                  subject_name, visit);
        end
    else
        error('Could not find Raw nirs data for %s - %s', subject_name, visit);
    end

    condition_aux = [visit '_origin_aux' ppl_tag];

    file_aux_raw = nst_get_bst_func_files(subject_name, condition_aux, 'Raw');
    if ~isempty(file_aux_raw)
        file_aux_raw = file_aux_raw{1};
        detect_heart_beats_ecg(file_aux_raw, options);
    else
        error('Could not find Raw aux data for %s - %s', subject_name, visit);
    end
end

end

function detect_heart_beats_ecg(sFile, options)
sEvents = in_bst_data(sFile, 'Events');
if ~ismember(options.heart_beats.event_name, {sEvents.Events.label})
    bst_process('CallProcess', 'process_evt_detect_multi_ecg', sFile, [], ...
                'channelnames', strjoin(options.heart_beats.ecg_channel_labels, ', '), ...
                'timewindow',   [], ...
                'dt',           0.05, ...
                'eventname_good_ecg',  options.heart_beats.good_ecg_event_name, ...
                'eventname',    options.heart_beats.event_name);
end
end

function detect_heart_beats_nirs(file_nirs, cardiac_channel_ref)
    data_tmp = in_bst_data(file_nirs , 'Events');
    if ~ismember('cardiac', {data_tmp.Events.label})
        % Process: Detect: cardiac
        bst_process('CallProcess', 'process_evt_nirs_cardiac', file_nirs, [], ...
                   'eventname',    'cardiac', ...
                    'channelname',  cardiac_channel_ref, ...
                   'timewindow',   [0, 875.04], ...
                   'bandpass',     [1, 7], ...
                   'threshold',    3, ...
                   'blanking',     0.12, ...
                   'isnoisecheck', 1);
    end
end
