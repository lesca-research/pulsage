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
        varargout{1} = imported_files;
        varargout{2} = redone;
        return;
    case 'load_markings'
        files_raw = {};
        files_fix = {};
        condition = ['origin' get_ppl_tag()];
        preproc_cond = ['preproc_' get_ppl_tag()];
        for iacq=1:length(arg1)
            acq = arg1(iacq);
            file_raw = nst_get_bst_func_files(acq.acq_tag, condition, 'Raw');
            if isempty(file_raw)
                error('Missing Raw data file for %s, %s\n', ...
                      acq.acq_tag, condition);
            end
            files_raw{end+1} = file_raw{1};

            file_fix = nst_get_bst_func_files(acq.acq_tag, preproc_cond, 'Raw_Chan_Fix');
            if isempty(file_fix)
                error('Missing Raw_Chan_Fix data file for %s, %s. Run preprocessing first.\n', ...
                      acq.acq_tag, preproc_cond);
            end 
            files_fix{end+1} = file_fix{1};

        end
        files_to_sync = [files_raw files_fix];
        
        create_dirs(options);
        sync_markings(files_to_sync, options, 'load_all');
        varargout{1} = files_to_sync;
        return;
    case 'save_markings'
        files_raw = {};
        files_fix = {};
        condition = ['origin' get_ppl_tag()];
        preproc_cond = ['preproc_' get_ppl_tag()];
        for iacq=1:length(arg1)
            acq = arg1(iacq);
            file_raw = nst_get_bst_func_files(acq.acq_tag, condition, 'Raw');
            if isempty(file_raw)
                error('Missing Raw data file for %s, %s\n', ...
                      acq.acq_tag, condition);
            end
            files_raw{end+1} = file_raw{1};

            file_fix = nst_get_bst_func_files(acq.acq_tag, preproc_cond, 'Raw_Chan_Fix');
            if isempty(file_fix)
                error('Missing Raw_Chan_Fix data file for %s, %s. Run preprocessing first.\n', ...
                      acq.acq_tag, preproc_cond);
            end 
            files_fix{end+1} = file_fix{1};

        end
        files_to_sync = [files_raw files_fix];
        
        create_dirs(options);
        sync_markings(files_to_sync, options, 'export_all');
        varargout{1} = files_to_sync;
        return;
    case 'preprocessing'
        preprocessing(arg1, options);
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
if options.pulsatility_az
    pulsatility_index_labels = {'Pulsatility_avg', 'Pulsatility_std', 'Pulsatility_min', 'Pulsatility_max'};
else
    pulsatility_index_labels = {'Pulsatility'};
end
ar_index_labels = [{'AR_Mx_mova', 'AR_Mx_block'} tfa_index_labels pulsatility_index_labels];
ar_indices = struct(); 
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


preproc_folder = sprintf('preproc_%s', get_ppl_tag());

all_features = dictionary();

for iacq=1:length(acq_defs)
    acq_name = acq_defs(iacq).acq_name;
    subject_tag = acq_defs(iacq).subject_tag;

    fprintf('Processing %s\n', acq_name);

    sFile_chan_fix = nst_get_bst_func_files(acq_name, preproc_folder, 'Raw_Chan_Fix');
    if isempty(sFile_chan_fix)
        warning('Cannot find "preproc/Raw_Chan_Fix" data for "%s".', acq_name);
        continue;
    elseif iscell(sFile_chan_fix) && length(sFile_chan_fix) > 1
        warning('Expect only 1 raw file "preproc/Raw_Chan_Fix", but %d found for "%s".',  ...
            length(sFile_chan_fix), acq_name);
        continue;
    end

    if iscell(sFile_chan_fix)
        sFile_chan_fix = sFile_chan_fix{1};
    end
        
    data_tmp = in_bst_data(sFile_chan_fix, 'Events', 'ChannelFlag');
    events = data_tmp.Events;

    [~, iStudy] = bst_get('AnyFile', sFile_chan_fix);
    sChannel = bst_get('ChannelForStudy', iStudy);
    channels = in_bst_channel(sChannel.FileName);
    channel_labels = {channels.Channel.Name};

    cfbi_idx_chans = find(ismember(channel_labels, CBFi_channel_labels) & (data_tmp.ChannelFlag==1)');
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

            spline_item = [acq_name '/' proc_folder '/' condition '_Spline' ];
            %TODO implement spline fit
            [sFile_spline, redone] = nst_run_bst_proc(spline_item, 0, ...
                                                    'process_psa_spline_fit', sFile_cond, [], ...
                                                    'channel_names', strjoin(CBFi_channel_labels, ', '), ...
                                                    'use_gap_events', 1, ...
                                                    'use_discard_events', 1, ...
                                                    'prior_interpolation_dt', options.spline_fit.prior_interpolation_dt, ...
                                                    'spline_knot_dt', options.spline_fit.knot_dt);
            sFiles_cond{end+1} = sFile_spline;
            prefixes{end+1} = [condition '_spline_'];
            suffixes{end+1} = ['_spline_' condition];
        end
    else
        raw_copy_item = [acq_name '/' proc_folder '/' 'Raw' ];
        [sFile_cond, redone] = nst_run_bst_proc(raw_copy_item, 0, ...
                                                'process_duplicate', sFile_chan_fix, [], ...
                                                'target', 1, ...  % Duplicate data files
                                                'tag',    '_copy');

        spline_item = [acq_name '/' proc_folder '/Spline' ];
        [sFile_spline, redone] = nst_run_bst_proc(spline_item, 0, ...
                                                    'process_psa_spline_fit', sFile_chan_fix, [], ...
                                                    'channel_names', CBFi_channel_labels, ...
                                                    'use_gap_events', 1, ...
                                                    'use_discard_events', 1, ...
                                                    'prior_interpolation_dt', options.spline_fit.prior_interpolation_dt, ...
                                                    'spline_knot_dt', options.spline_fit.knot_dt);

        sFiles_cond = {sFile_cond, sFile_spline};
        prefixes = {'' , 'spline_'};
        suffixes = {'', '_spline'};
    end

    for ifile=1:length(sFiles_cond)
        
        prefix = prefixes{ifile};
        suffix = suffixes{ifile};

        fgap_item = [acq_name '/' proc_folder '/' prefix 'gap_filled'];
        [sFile_fgaps, redone] = nst_run_bst_proc(fgap_item, 0, ...
                                                'process_psa_fill_gaps', sFiles_cond{ifile}, [], ...
                                                'option_channels',    [CBFi_channel_labels {BP_channel_label}], ...
                                                'option_ar_win_size', options.fill_gaps.ar_win_size_sec, ...
                                                'option_ar_order',    options.fill_gaps.ar_order_sec);
    
        % TODO email laurette pour stim cog
        curated_AR_item = [acq_name '/' proc_folder '/' prefix 'curated_for_AR' ];
        discard_BP_event = options.BP_channel_label;
        [sFile_curated_AR, redone] = nst_run_bst_proc(curated_AR_item, 0, ...
                                                    'process_psa_curate_event', sFile_fgaps, [], ...
                                                    'use_discard_channel_events',   1, ...
                                                    'discard_global_events', [discard_BP_event ', extrasystole'], ...
                                                    'channelnames',  strjoin([channel_labels(cfbi_idx_chans) {BP_channel_label}], ','));        
        curated_puls_item = [acq_name '/' proc_folder '/' prefix 'curated_for_puls' ];
        [sFile_curated_puls, redone] = nst_run_bst_proc(curated_puls_item, 0, ...
                                        'process_psa_curate_event', sFile_fgaps, [], ...
                                        'use_discard_channel_events',   1, ...
                                        'discard_global_events', 'discard_cardiac, extrasystole', ...
                                        'channelnames',  strjoin(channel_labels(cfbi_idx_chans), ','));        

        if options.pulsatility.do
            pulsatility_item = [acq_name '/' proc_folder '/' prefix 'pulsatility' ];
            
            [sFile_puls, redone] = nst_run_bst_proc(pulsatility_item, 0, ...
                                                    'process_psa_pulsatility', sFile_curated_puls, [], ...
                                                    'timewindow',    [], ...
                                                    'channelnames',  strjoin(channel_labels(cfbi_idx_chans), ','), ...
                                                    'heart_beat_event', options.heart_beats.event_name, ...
                                                    'function', 'pulse_amplitude_normed'); 

            pulsatility_avg_item = [acq_name '/' proc_folder '/' prefix 'pulsatility_avg' ];
            [sFile_puls_avg, redone] = nst_run_bst_proc(pulsatility_avg_item, 0, ...
                                                        'process_psa_average_time_robust', sFile_puls, [], ...
                                                        'pct_filter_range', [5, 95], ...
                                                        'use_discard_events',   0, ...
                                                        'avg_func',  'mean', ...  % Arithmetic average:  mean(x)
                                                        'overwrite',        0);

            pulsatility_avg_data = in_bst(sFile_puls_avg);

            for chan_idx=cfbi_idx_chans
                if pulsatility_avg_data.ChannelFlag(chan_idx) == 1
                    pi = pulsatility_avg_data.F(chan_idx, 1);
                    tcd_chan_label = channel_labels{chan_idx};
                    pi_label = protect_field_label(['Pulsatility_' tcd_chan_label suffix]);
                    key = {{subject_tag, pi_label}};
                    if isConfigured(all_features) && isKey(all_features, key)
                        error('%s already computed for %s, duplicate condition or duplicate entry in manifest?', pi_label, subject_tag);
                    end
                    all_features(key) = pi;
                end
            end
        end

        bp_channel_idx = ismember(channel_labels, BP_channel_label);
        dataMat = in_bst_data(sFile_curated_AR, 'ChannelFlag', 'Events', 'Time');

        if all(bp_channel_idx==0) || dataMat.ChannelFlag(bp_channel_idx) == -1
            warning('No BP channel available in %s, autoregulation not computed.', ...
                    acq_name)
            continue;
        end

        events = dataMat.Events;
        idiscard_evt = find(strcmp('discard_BP', {events.label}));
        if ~isempty(idiscard_evt) && sum(diff(events(idiscard_evt).times)) / diff(dataMat.Time([1 end])) > 0.5
            warning('Ignore BP %s signal because more than 50%% is discarded', ...
                    condition);
            continue;
        end

        mx_item = [acq_name '/' proc_folder '/' prefix 'Mx' ];
        % TODO expose otions for pct filter
        [sFile_mx, redone] = nst_run_bst_proc(mx_item, 0, ...
                                              'process_psa_mx', sFile_curated_AR, [], ...
                                              'abp_channel', BP_channel_label, ...
                                              'cbf_channels', channel_labels(cfbi_idx_chans));
        if  options.low_pass_filter.do
            % Low-pass filter
            filter_item = [acq_name '/' proc_folder '/' prefix 'filtered' ];
            [sFile_filtered, redone] = nst_run_bst_proc(filter_item, 0, ...
                'process_bandpass', sFile_curated_AR, [], ...
                'sensortypes', 'BP, TCD', ...
                'highpass',    0, ...
                'lowpass',     options.low_pass_filter.low_cutoff, ...
                'tranband',    0, ...
                'attenuation', 'strict', ...  % 60dB
                'ver',         '2019', ...  % 2019
                'mirror',      0, ...
                'overwrite',   0);
        else
            sFile_filtered = sFile_curated_AR;
        end

        if options.moving_average.do
            win_avg_comment = [prefix 'mov_average'];
            win_avg_item = [acq_name '/' proc_folder '/' win_avg_comment ];

            [sFile_mov_avg, redone] = nst_run_bst_proc(win_avg_item, 0, ...
                                                        'process_psa_moving_average', sFile_filtered, [], ...
                                                        'option_channels', [BP_channel_label channel_labels(cfbi_idx_chans)], ...
                                                        'option_win_size', options.moving_average_window_sec);
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
        else
            sFile_mov_avg = sFile_filtered;
        end

        do_tfa = 0;
        if do_tfa
            tfa_item = [acq_name '/' proc_folder '/' prefix 'TFA' ];
            [sFile_tfa, redone] = nst_run_bst_proc(tfa_item, 0, ...
                                                    'process_psa_tfa', sFile_mov_avg, [], ...
                                                    'ref_channel', BP_channel_label, ...
                                                    'transfer_channels', channel_labels(cfbi_idx_chans));
            tfa_data = in_bst(sFile_tfa);
        end
        MatFlag = in_bst_data(sFile_fgaps, 'ChannelFlag');

        mx_data = in_bst(sFile_mx);
        for i_tcd_chan=1:length(cfbi_idx_chans)
            tcd_chan_label = channel_labels{cfbi_idx_chans(i_tcd_chan)};

            mx_index_label = protect_field_label(['AR_Mx_block_' tcd_chan_label]);
            i_mx = strcmp(mx_index_label, mx_data.Description);
            ar_index_label = [mx_index_label suffix];
            key = {{subject_tag, ar_index_label}};

            if MatFlag.ChannelFlag(cfbi_idx_chans(i_tcd_chan)) == 1
                if isConfigured(all_features) &&  isKey(all_features, key)
                    error('%s already computed for %s, duplicate condition or duplicate entry in manifest?', ar_index_label, subject_tag);
                end
                all_features(key) = mx_data.Value(i_mx);
            end
        end
    end % end of loop over conditions (eg rest, stand)
end % end of loop over acquisitions

feature_table = dict_to_table(all_features, {'subject_id'});
outputFile = fullfile(options.result_dir, 'hemodynamics.xlsx');
writetable(feature_table, outputFile, 'Sheet', 'Hemodynamics_Features');
disp(feature_table);
fprintf('Wrote %d rows and %d columns to %s\n', ...
        height(feature_table), width(feature_table), outputFile);
end

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

function preprocessing(acq_defs, options)

for iacq=1:length(acq_defs)
    acq_name = acq_defs(iacq).acq_name;

    fprintf('Preprocessing %s\n', acq_name);

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

    add_default_events(bst_process('GetInputStruct', file_raw), options);

    preproc_folder = sprintf('preproc_%s', get_ppl_tag());

    scale_fix_item = [acq_name '/' preproc_folder '/Raw_Scale_Fix' ];
    if isfield(acq_defs, 'TCD_scaling') && ~isempty(acq_defs(iacq).TCD_scaling)
        [sFile_scale_fix, redone] = nst_run_bst_proc(scale_fix_item, 0, ...
                                                    'process_scale', file_raw, [], ...
                                                    'factor',      acq_defs(iacq).TCD_scaling, ...
                                                    'sensortypes', 'TCD', ...
                                                    'overwrite',   0);
    else
        [sFile_scale_fix, redone] = nst_run_bst_proc(scale_fix_item, 0, ...
                                                    'process_duplicate', file_raw, [], ...
                                                    'target', 1, ...  % Duplicate data files
                                                    'tag',    '_copy');
    end

    % Add extra channels for spline-fitting preview
    [channelsMat, iStudy] = load_channels(sFile_scale_fix);
    DataMat = in_bst_data(sFile_scale_fix);
    cbfi_idx_chans = find(ismember({channelsMat.Channel.Name}, options.CBFi_channel_labels));
    if ~ismember([options.CBFi_channel_labels{1} spline_suffix()], {channelsMat.Channel.Name})
        ichan_denoised = length(channelsMat.Channel) + 1;
        DataMat.ChannelFlag = [DataMat.ChannelFlag ; zeros(length(cbfi_idx_chans), 1)];
        for cbfi_chan_idx=cbfi_idx_chans
            channel_def = channelsMat.Channel(cbfi_chan_idx);
            channel_def.Name = [channel_def.Name spline_suffix()];
            channelsMat.Channel = [channelsMat.Channel channel_def];
            DataMat.ChannelFlag(ichan_denoised) = DataMat.ChannelFlag(cbfi_chan_idx);
            DataMat.F = [DataMat.F ; zeros(1, size(DataMat.F, 2)) + ichan_denoised];
            ichan_denoised = ichan_denoised + 1;
        end
        db_set_channel(iStudy, channelsMat, 2, 0);
        bst_save(file_fullpath(sFile_scale_fix), DataMat, 'v7');
    end

    chan_fix_item = [acq_name '/' preproc_folder '/Raw_Chan_Fix' ];
    if isfield(acq_defs, 'channels_to_swap') && ~isempty(acq_defs(iacq).channels_to_swap)
        [sFile_chan_fix, redone] = nst_run_bst_proc(chan_fix_item, 0, ...
                                                    'process_psa_swap_channels', sFile_scale_fix, [], ...
                                                    'channelname_1',   acq_defs(iacq).channels_to_swap{1}, ...
                                                    'channelname_2',   acq_defs(iacq).channels_to_swap{2});
    else
        [sFile_chan_fix, redone] = nst_run_bst_proc(chan_fix_item, 0, ...
                                                    'process_duplicate', sFile_scale_fix, [], ...
                                                    'target', 1, ...  % Duplicate data files
                                                    'tag',    '_copy2');
    end

    DataMat = in_bst_data(sFile_chan_fix);
    if size(DataMat.ChannelFlag, 1) ~= length(channelsMat.Channel)
        ichan_denoised = length(DataMat.ChannelFlag) + 1;
        DataMat.ChannelFlag = [DataMat.ChannelFlag ; zeros(length(cbfi_idx_chans), 1)];
        for cbfi_chan_idx=cbfi_idx_chans
            DataMat.ChannelFlag(ichan_denoised) = DataMat.ChannelFlag(cbfi_chan_idx);
            DataMat.F = [DataMat.F ; zeros(1, size(DataMat.F, 2)) + ichan_denoised];
            ichan_denoised = ichan_denoised + 1;
        end
        bst_save(file_fullpath(sFile_chan_fix), DataMat, 'v7');
    end
    
    % TODO crop here to avoid too long signals for spline fitting

    dest_channels = cellfun(@(s) [s spline_suffix()], options.CBFi_channel_labels, 'UniformOutput', false);
    bst_process('CallProcess', 'process_psa_spline_fit', sFile_chan_fix, [], ...
                'channel_names', strjoin(options.CBFi_channel_labels, ', '), ...
                'spline_channel_names', strjoin(dest_channels, ', '), ...
                'ignore_gap_events', 1, ...
                'ignore_discard_events', 1, ...
                'prior_interpolation_dt', options.spline_fit.prior_interpolation_dt, ...
                'spline_knot_dt', options.spline_fit.knot_dt);

    detect_heart_beats(sFile_chan_fix, ...
                       options.heart_beats.ecg_channel_labels, ...
                       options.heart_beats.event_name);

    gap_bp_label = ['gap_' options.BP_channel_label];
    detect_BP_gaps(sFile_chan_fix, options.BP_channel_label, gap_bp_label);

    add_default_events(bst_process('GetInputStruct', sFile_chan_fix), options);
end

end

function suffix = spline_suffix()
suffix = '_denoised';
end

function detect_BP_gaps(sFile, BP_channel_label, BP_gap_event_label)
sEvents = in_bst_data(sFile, 'Events');
evt_gap = ismember({sEvents.Events.label}, BP_gap_event_label);

[~, iStudy] = bst_get('AnyFile', sFile);
sChannel = bst_get('ChannelForStudy', iStudy);
channels = in_bst_channel(sChannel.FileName);
channel_names = {channels.Channel.Name};

if ismember(BP_channel_label, channel_names) && (~any(evt_gap) || isempty(sEvents.Events(evt_gap).times))
    bst_process('CallProcess', 'process_evt_detect_artifacts', sFile, [], ...
                'channelname', BP_channel_label, ...
                'signal_type',  'BP', ...
                'event_name',    BP_gap_event_label);
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

options.discard_segments = struct([]);

options.fill_gaps.do = 1;
options.fill_gaps.ar_win_size_sec = 15;
options.fill_gaps.ar_order_sec = 1;

options.spline_fit.prior_interpolation_dt = 0.03; % sec
options.spline_fit.knot_dt = 0.15; % sec

options.low_pass_filter.do = 0;
options.low_pass_filter.low_cutoff = 0.2; % Hz

options.moving_average.do = 0;
options.moving_average_window_sec = 10;

options.pulsatility.do = 1;
options.pulsatility_az = 0;

options.deglitch.do = 0;
options.deglitch.redo = 0;
options.deglitch.agrad_std_factor = 2.5;

%fig OPTIONS
options.make_figs = 1; 
options.save_fig_method = 'saveas'; % 'saveas', 'export_fig'
options.export_fig_dpi = 90;
options.fig_dir = ''; % the directory where 'pulsatility_figures' will be generated (default in data_analysis)
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
    add_default_events(sFiles, options);
    sFiles = bst_process('CallProcess', 'process_set_comments', sFiles, [], ...
                         'tag', 'Raw', 'isindex', 0);
    files_in = [files_in sFiles.FileName];
end

end

function add_default_events(sFile, options)
warning off MATLAB:load:variableNotFound
data_tmp = in_bst_data(sFile.FileName, 'Events');
warning on MATLAB:load:variableNotFound
if ~isempty(data_tmp.Events)
    data_event_labels = {data_tmp.Events.label};
else
    data_event_labels = {};
end

[~, iStudy] = bst_get('AnyFile', sFile.FileName);
sChannel = bst_get('ChannelForStudy', iStudy);
channels = in_bst_channel(sChannel.FileName);
channel_names = {channels.Channel.Name};

%TODO save all markings again

cbfi_idx_chans = find(ismember(channel_names, options.CBFi_channel_labels));
BP_chan_idx = find(ismember(channel_names, options.BP_channel_label));

tag_events = db_template('event');
ievt = 1;
for chan_idx=[cbfi_idx_chans BP_chan_idx]
    gap_event_label = ['gap_' channel_names{chan_idx}];
    if ~ismember(gap_event_label, data_event_labels)
        tag_events(ievt).label = gap_event_label;
        tag_events(ievt).times = zeros(2,0);
        ievt = ievt + 1;
    end

    discard_event_label = ['discard_' channel_names{chan_idx}];
    if ~ismember(discard_event_label, data_event_labels)
        tag_events(ievt).label = discard_event_label;
        tag_events(ievt).times = zeros(2,0);
        ievt = ievt + 1;
    end
end

if ~ismember('extrasystole', data_event_labels)
    tag_events(ievt).label = 'extrasystole';
    tag_events(ievt).times = zeros(2,0);
    ievt = ievt + 1;
end

if ~ismember('discard_cardiac', data_event_labels)
    tag_events(ievt).label = 'discard_cardiac';
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
    process_nst_import_csv_events('import_events', [], sFile, tag_events);
end
end

function create_dirs(options)
create_dir(options.result_dir);
create_dir(options.fig_dir);
create_dir(options.export_dir_events);
create_dir(options.export_dir_channel_flags);
create_dir(options.export_dir_preprocessed);
end

function sync_markings(sFiles, options, sync_mode)

gap_event_labels = [{['gap_' options.BP_channel_label]} ...
                    cellfun(@(s) ['gap_' s], options.CBFi_channel_labels, 'UniformOutput', false)];
discard_event_labels = [{['discard_' options.BP_channel_label]} ...
                        cellfun(@(s) ['discard_' s], options.CBFi_channel_labels, 'UniformOutput', false)];
io_options.events_to_sync = [options.conditions ...
                             {'extrasystole', 'discard_cardiac', options.heart_beats.event_name} ...
                             gap_event_labels discard_event_labels];

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
        add_default_events(sInputs(ifile), options);
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

function [channelMat, iStudy] = load_channels(sFile)
[~, iStudy] = bst_get('AnyFile', sFile);
sChannel = bst_get('ChannelForStudy', iStudy);
channelMat = in_bst_channel(sChannel.FileName);
end
