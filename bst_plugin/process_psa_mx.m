function varargout = process_psa_mx( varargin )

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
sProcess.Comment     = 'Mean Flow Index';
sProcess.FileTag     = 'Mx';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Frequency';
sProcess.Index       = 3; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'matrix', 'matrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options
sProcess.options.abp_channel.Comment = 'ABP Channel: ';
sProcess.options.abp_channel.Type    = 'text';
sProcess.options.abp_channel.Value   = '';

sProcess.options.cbf_channels.Comment = 'CBF-related Channels (comma-separated types or names): ';
sProcess.options.cbf_channels.Type    = 'text';
sProcess.options.cbf_channels.Value   = '';

% Minimum delay between events
sProcess.options.block_duration.Comment = 'Block duration: ';
sProcess.options.block_duration.Type    = 'value';
sProcess.options.block_duration.Value   = {10, 'sec', 0};

sProcess.options.pct_filter.Comment = 'Percentage of extreme values to filter: ';
sProcess.options.pct_filter.Type    = 'value';
sProcess.options.pct_filter.Value   = {5, '%', 0};

sProcess.options.max_nan_pct.Comment = 'Percentage of NaN values allowed per block: ';
sProcess.options.max_nan_pct.Type    = 'value';
sProcess.options.max_nan_pct.Value   = {25, '%', 0};

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
    
    channels = in_bst_channel(sInputs(iInput).ChannelFile);
        
    abp_sig_idx_chan = channel_find(channels.Channel, sProcess.options.abp_channel.Value);
    cbf_sig_idx_chan = channel_find(channels.Channel, sProcess.options.cbf_channels.Value);
    
    ref_signal = sDataIn.F(abp_sig_idx_chan, :);
    dt = diff(sDataIn.Time(1:2));
    for ichan=1:length(cbf_sig_idx_chan)
        idx_chan = cbf_sig_idx_chan(ichan);
        cbf_signal = sDataIn.F(idx_chan, :);      
        block_size = round(sProcess.options.block_duration.Value{1} / dt);
        mx_out = Compute(ref_signal, cbf_signal, block_size, ...
                         sProcess.options.pct_filter.Value{1}/100, ...
                         sProcess.options.max_nan_pct.Value{1}/100);
        chan_name = channels.Channel(idx_chan).Name;
        mx_values.(protect_field_label(['AR_Mx_block_' chan_name])) = mx_out;
    end   
    
    % Save as matrix
    comment = [sInputs(iInput).Comment  ' | Mx' ];
    OutputFiles{iInput} = nst_save_table_in_bst(struct2table(mx_values), sInputs(iInput).SubjectName, sInputs(iInput).Condition,...
                                                 comment);
end

end

function label = protect_field_label(label)
label = strrep(label, '.', '_dot_');
end


%% ===== Compute =====
function [Mx] = Compute(bp_signal, cbfi_signal, block_size, pct_filter, max_nan_prop)

if size(bp_signal, 1) > size(bp_signal, 2)
    block_window = [block_size, 1];
else
    block_window = [1, block_size];
end

% fprintf('Compute blocks\n');
bp_signal_blockproc = blockproc(bp_signal, block_window, @(x)mean(x.data));
bp_signal_block_test = block_avg(bp_signal, block_size, -1);
assert(isequaln(bp_signal_block_test, bp_signal_blockproc));

bp_signal_block = block_avg(bp_signal, block_size, 1-max_nan_prop);
cbfi_signal_block = block_avg(cbfi_signal, block_size, 1-max_nan_prop);

m_nan = isnan(bp_signal_block) | isnan(cbfi_signal_block);
bp_signal_block = bp_signal_block(~m_nan);
cbfi_signal_block = cbfi_signal_block(~m_nan);

% fprintf('Filter blocks\n');
[bp_signal_block_filtered, cbfi_signal_block_filtered] = pct_filter2(bp_signal_block, cbfi_signal_block, pct_filter);

% TODO expose minimum number of blocks
if length(bp_signal_block_filtered) >= 10
    % fprintf('Compute block correlations\n');
    [r, p] = corrcoef(bp_signal_block_filtered, cbfi_signal_block_filtered);
    
    if p > 0.05
        warning('No significant correlation between block signals of BP and CBFi (r=%1.2f)', r(1,2));
    end
    
    Mx = r(1,2);
else
    warning('No enough blocks (%d) to compute Mx', length(bp_signal_block_filtered));
    Mx = NaN;
end
end

function block_signal = block_avg(signal, block_size, prop_min_valid_values)
nb_samples = length(signal);
nb_blocks = ceil(nb_samples / block_size);
block_signal = nan(1, nb_blocks);
for iblock=1:nb_blocks
    block_smpl_start = ((iblock-1) * block_size ) + 1;
    block_smpl_end = min(nb_samples, block_smpl_start + block_size - 1);
    chunk = signal(block_smpl_start:block_smpl_end);
    if prop_min_valid_values == -1
        block_signal(iblock) = mean(chunk);
    elseif sum(~isnan(chunk))/block_size >= prop_min_valid_values
        if any(isnan(chunk))
            fprintf('%d nans ignored in block %d\n', sum(isnan(chunk)), iblock);
        end
        block_signal(iblock) = mean(chunk, "omitnan");
    end
end
end

function [sig1_filtered, sig2_filtered] = pct_filter2(sig1, sig2, pct)

[~, sig1_sorted_idx] = sort(sig1);
sig1_low_idx = max(1, round(length(sig1)*pct/2));
sig1_high_idx = length(sig1) - sig1_low_idx;

sig1 = sig1(sig1_sorted_idx(sig1_low_idx:sig1_high_idx));
sig2 = sig2(sig1_sorted_idx(sig1_low_idx:sig1_high_idx));

[~, sig2_sorted_idx] = sort(sig2);
sig2_low_idx = max(1, round(length(sig2)*pct/2));
sig2_high_idx = length(sig2) - sig2_low_idx;

sig1_filtered = sig1(sig2_sorted_idx(sig2_low_idx:sig2_high_idx));
sig2_filtered = sig2(sig2_sorted_idx(sig2_low_idx:sig2_high_idx));

end


function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end
