function varargout = process_grouped_analysis( varargin )
% GROUP correlation: TCD PI (TSV) vs NIRS channel PI (Brainstorm inputs)
%
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill
% University
% This software is distributed under the terms of the GNU General Public
% License (GPLv3): http://www.gnu.org/copyleft/gpl.html
% FOR RESEARCH PURPOSES ONLY - PROVIDED "AS IS"
% =============================================================================@
%
% Authors: Aymen Zire, 2025
%
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'GROUP: TCD PI ↔ NIRS PI (corr)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1651;
    sProcess.Description = 'Loads /home/lesca-student/PULSATILITY_TCD.tsv and correlates with NIRS channels.';

    % Inputs: one file per subject (must contain NIRS PI per channel)
    sProcess.InputTypes  = {'data'};      % or 'pdata'
    sProcess.OutputTypes = {'pdata'};     % stat-at-sensors
    sProcess.nMinFiles   = 2;
    sProcess.nInputs     = 1;

    % === Options
    sProcess.options.tsv_path.Comment = 'TSV path:';
    sProcess.options.tsv_path.Type    = 'text';
    sProcess.options.tsv_path.Value   = '/home/lesca-student/PULSATILITY_TCD.tsv';

    % Which TCD column to use
    sProcess.options.tcd_column.Comment = { ...
        'R_rest_avg','R_stand_avg','L_rest_avg','L_stand_avg','R','L',''; ...
        'Pulsatility_TCD_MCA_R_rest_avg','Pulsatility_TCD_MCA_R_stand_avg', ...
        'Pulsatility_TCD_MCA_L_rest_avg','Pulsatility_TCD_MCA_L_stand_avg', ...
        'R','L',''};
    sProcess.options.tcd_column.Type    = 'radio_linelabel';
    sProcess.options.tcd_column.Value   = 'Pulsatility_TCD_MCA_R_rest_avg';

    % Correlation method
    sProcess.options.method.Comment  = {'Pearson','Spearman',''; 'pearson','spearman',''};
    sProcess.options.method.Type     = 'radio_linelabel';
    sProcess.options.method.Value    = 'pearson';

    % Missing values policy
    sProcess.options.nan_policy.Comment = {'Pairwise (omit NaN)','Complete cases only',''; ...
                                           'pairwise','complete',''};
    sProcess.options.nan_policy.Type  = 'radio_linelabel';
    sProcess.options.nan_policy.Value = 'pairwise';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    %Options
    tsv_path   = sProcess.options.tsv_path.Value;
    tcd_colkey = sProcess.options.tcd_column.Value;
    method     = sProcess.options.method.Value;
    nan_policy = sProcess.options.nan_policy.Value;

    % Read TSV
    assert(exist(tsv_path,'file')==2, ['TSV not found: ' tsv_path]);
    t = readtable(tsv_path,'FileType','text','Delimiter','\t','ReadVariableNames',true);
    assert(any(strcmpi(t.Properties.VariableNames,'Subject')), 'TSV must contain column "Subject"');

    % Map chosen key -> actual column
    tcd_colname = map_tcd_colname(t,tcd_colkey);
    assert(any(strcmp(t.Properties.VariableNames, tcd_colname)), ['Column not found in TSV: ' tcd_colname]);

    subj_tsv = string(t.Subject);
    vals_tcd = t.(tcd_colname);

    % =subject IDs from Brainstorm inputs 
    nIn = numel(sInputs);
    subj_in = strings(nIn,1);

    for i = 1:nIn
        subj = sInputs(i).SubjectName;
        if isempty(subj)
            [sSubj, ~] = bst_get('SubjectForStudy', sInputs(i).iStudy);
            if ~isempty(sSubj) && isfield(sSubj,'Name') && ~isempty(sSubj.Name)
                subj = sSubj.Name;
            else
                error('Cannot get subject name for: %s', sInputs(i).FileName);
            end
        end
        subj_in(i) = string(subj);
    end

    % Truncate to "..._T0" so it matches the TSV "Subject" key
    for i = 1:nIn
        if contains(subj_in(i), '_T0')
            subj_in(i) = extractBefore(subj_in(i), '_T0') + "_T0";
        end
    end

    % Align inputs ↔ TSV rows
    idx_tsv = nan(nIn,1);
    for i=1:nIn
        j = find(subj_tsv == subj_in(i), 1, 'first');
        if ~isempty(j)
            idx_tsv(i) = j;
        end
    end
    keep_inputs = ~isnan(idx_tsv);
    nMatched    = sum(keep_inputs);
    if nMatched < 2
        error('Not enough matched subjects between inputs and TSV to compute correlations. Matched=%d', nMatched);
    end

    %    common channels + stack NIRS PI
    [ChannelMat, ChannelList] = getCommonChannels(sInputs);
    nChan    = numel(ChannelList);
    all_nirs = nan(nIn, nChan);

    for i=1:nIn
        data_i = getDataForChannels(sInputs(i), ChannelList);
        assert(numel(data_i.F)==nChan, 'Mismatch channels vs data.F for input %d', i);
        all_nirs(i,:) = data_i.F(:).';
    end

    % Keep matched rows only
    all_nirs = all_nirs(keep_inputs,:);
    tcd_vec  = vals_tcd(idx_tsv(keep_inputs));

    %Channel-wise correlation
    [r_map, p_map, n_map] = corr_per_channel(all_nirs, tcd_vec, method, nan_policy);

    % -Save as Brainstorm pdata (stat at sensors)
    sOutput = db_template('statmat');
    sOutput.Type         = 'pdata';
    sOutput.Time         = 1;                          % one dummy time sample (UI expects an index)
    sOutput.nComponents  = 1;
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 'r';
    sOutput.ChannelFlag  = ones(1,nChan);
    sOutput.GoodChannel  = ones(1,nChan);              % helps sensors display
    sOutput.Options.SensorTypes = 'NIRS';
    sOutput.Method       = ['channelwise-' method];
    sOutput.Measure      = 'corr';
    sOutput.RowNames     = ChannelList(:);
    sOutput.Comment      = sprintf('Corr(TCD:%s, NIRS PI) [%s, %s]', tcd_colkey, method, nan_policy);

    % Store stats
    sOutput.tmap = r_map(:);                           % r per channel
    sOutput.pmap = p_map(:);                           % p per channel
    sOutput.df   = max(n_map-2, 0)';                   % df ~ n-2

    % Save under Group analysis / TCD_NIRS
    out_subject   = 'Group analysis';
    out_condition = 'TCD_NIRS';
    [sSubject, iSubject] = bst_get('Subject', out_subject, 1);
    if isempty(sSubject), [sSubject, iSubject] = db_add_subject(out_subject, [], 1, 0); end
    iStudy = db_add_condition(out_subject, out_condition);

    % Ensure a matching channel file in that study (same channels/order)
    [~, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
    db_set_channel(iChannelStudy, ChannelMat, 0, 0);

    % Sanity: channel file must match vector sizes
    sStudy  = bst_get('Study', iStudy);
    chanMat = in_bst_channel(sStudy.Channel.FileName);
    assert(numel(chanMat.Channel) == nChan, ...
        'Channel file size (%d) mismatch with tmap/pmap size (%d).', numel(chanMat.Channel), nChan);

    % Save
    outFile = bst_process('GetNewFilename', fileparts(sStudy.FileName), 'pdata_tcd_nirs_corr');
    bst_save(outFile, sOutput, 'v6');
    db_add_data(iStudy, outFile, sOutput);

    % Return to pipeline
    OutputFiles = {outFile};
end

%% ===== Helpers =====
function tcd_colname = map_tcd_colname(t, key)
    switch key
        case 'Pulsatility_TCD_MCA_R_rest_avg'
            tcd_colname = 'Pulsatility_TCD_MCA_R_rest_avg';
        case 'Pulsatility_TCD_MCA_R_stand_avg'
            tcd_colname = 'Pulsatility_TCD_MCA_R_stand_avg';
        case 'Pulsatility_TCD_MCA_L_rest_avg'
            tcd_colname = 'Pulsatility_TCD_MCA_L_rest_avg';
        case 'Pulsatility_TCD_MCA_L_stand_avg'
            tcd_colname = 'Pulsatility_TCD_MCA_L_stand_avg';
        case 'R'
            tcd_colname = 'R';
        case 'L'
            tcd_colname = 'L';
        otherwise
            error('Unknown key for TCD column: %s', key);
    end
end

function [r_map, p_map, n_map] = corr_per_channel(all_nirs, tcd_vec, method, nan_policy)
    nChan = size(all_nirs,2);
    r_map = nan(1,nChan);
    p_map = nan(1,nChan);
    n_map = nan(1,nChan);
    for c = 1:nChan
        x = all_nirs(:,c);
        y = tcd_vec(:);

        switch nan_policy
            case 'pairwise'
                ok = ~(isnan(x) | isnan(y));
            case 'complete'
                ok = ~(isnan(x) | isnan(y));   % same effect here (2 vars only)
            otherwise
                ok = ~(isnan(x) | isnan(y));
        end
        x = x(ok); y = y(ok);

        if numel(x) >= 3
            switch method
                case 'pearson'
                    [r,p] = corr(x, y, 'Type','Pearson', 'Rows','complete');
                case 'spearman'
                    [r,p] = corr(x, y, 'Type','Spearman', 'Rows','complete');
                otherwise
                    error('Unknown correlation method: %s', method);
            end
            r_map(c) = r;
            p_map(c) = p;
            n_map(c) = numel(x);
        else
            r_map(c) = NaN;
            p_map(c) = NaN;
            n_map(c) = numel(x);
        end
    end
end

function [ChannelMat,ChannelList]=getCommonChannels(sInputs)
    nb_inputs=length(sInputs);
    ChannelMat1 = in_bst_channel(sInputs(1).ChannelFile);
    ChannelList = getChannelsNames(ChannelMat1);
    for iInput=2:nb_inputs
        ChannelMat_i = in_bst_channel(sInputs(iInput).ChannelFile);
        tmpChannelList = getChannelsNames(ChannelMat_i);
        ChannelList = intersect(ChannelList, tmpChannelList, 'stable');
    end
    % rebuild ChannelMat with only common channels (preserve order)
    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    new_channel = repmat(db_template('channeldesc'), 1, 0);
    for i_old=1:length(ChannelMat.Channel)
        if ismember(ChannelMat.Channel(i_old).Name, ChannelList)
            new_channel(end+1) = ChannelMat.Channel(i_old); %#ok<AGROW>
        end
    end
    ChannelMat.Channel = new_channel;
    ChannelMat.Comment = ['NIRS sensors (' int2str(length(ChannelList)) ')'];
end

function ChannelsNames=getChannelsNames(channelMat)
    [nirs_ichans, ~] = channel_find(channelMat.Channel, 'NIRS');
    ChannelsNames = arrayfun(@(i) channelMat.Channel(i).Name, nirs_ichans, 'UniformOutput', false);
end

function data=getDataForChannels(sInput,ChannelList)
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    names_all = getChannelsNames(ChannelMat);
    data = in_bst_data(sInput.FileName);

    F = nan(numel(ChannelList),1);
    for i=1:numel(names_all)
        [C,~,ib] = intersect(names_all(i), ChannelList);
        if ~isempty(C)
            F(ib) = data.F(i);  
        end
    end
    data.F = F;
    data.ChannelFlag = ones(numel(ChannelList),1);
end
