function varargout = process_evt_detect_artifacts( varargin )
% PROCESS_EVT_DETECT_ARTIFACTS: Detect artifacts"

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
% Authors: Francois Tadel, Elizabeth Bock, 2011-2016, 
% Thomas Vincent 2025-2026

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Detect artifacts';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 44;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;


    sProcess.options.channelname.Comment = 'Channel name:';
    sProcess.options.channelname.Type    = 'text';
    sProcess.options.channelname.Value   = '';

    % Event prefix
    sProcess.options.event_name.Comment = 'Event name: ';
    sProcess.options.event_name.Type    = 'text';
    sProcess.options.event_name.Value   = 'artifacts';

    sProcess.options.signal_type.Comment = {'ECG', ...
                                            'Blood Pressure'; ...
                                            'ECG', 'BP'};
    sProcess.options.signal_type.Type    = 'radio_label';
    sProcess.options.signal_type.Value   = 'ECG';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>   

for iInput=1:length(sInputs)
    % Load the raw file descriptor
    isRaw = strcmpi(sInputs(iInput).FileType, 'raw');
    if isRaw
        DataMat = in_bst_data(sInputs(iInput).FileName, 'F', 'Time');
        sFile = DataMat.F;
    else
        DataMat = in_bst_data(sInputs(iInput).FileName, 'Time');
        sFile = in_fopen(sInputs(iInput).FileName, 'BST-DATA');
    end
    % Load channel file
    ChannelMat = in_bst_channel(sInputs(iInput).ChannelFile);

    if ~isempty(sProcess.options.channelname.Value)
        idx_chan = channel_find(ChannelMat.Channel, sProcess.options.channelname.Value);
        if isempty(idx_chan)
            error('Channel %s not found', sProcess.options.channelname.Value);
        end
    else
        error('No channel specified');
    end

    [F, TimeVector] = in_fread(sFile, ChannelMat, 1, [], idx_chan);
    artifact_times = Compute(F, TimeVector, sProcess.options.signal_type.Value);

    % ===== CREATE EVENTS =====
    % Basic events structure
    if ~isfield(sFile, 'events') || isempty(sFile.events)
        sFile.events = repmat(db_template('event'), 0);
    end
    
    % Event name
    newName = sProcess.options.event_name.Value;
    % Get the event to create
    iEvt = find(strcmpi({sFile.events.label}, newName));
    % Existing event: reset it
    if ~isempty(iEvt)
        sEvent = sFile.events(iEvt);
        sEvent.epochs     = [];
        sEvent.times      = [];
        sEvent.reactTimes = [];
    % Else: create new event
    else
        % Initialize new event
        iEvt = length(sFile.events) + 1;
        sEvent = db_template('event');
        sEvent.label = newName;
        % Get the default color for this new event
        sEvent.color = panel_record('GetNewEventColor', iEvt, sFile.events);
    end
    % Times, samples, epochs
    sEvent.times    = artifact_times;
    sEvent.epochs   = ones(1, size(sEvent.times,2));
    sEvent.channels = [];
    sEvent.notes    = [];
    % Add to events structure
    sFile.events(iEvt) = sEvent;

    % ===== SAVE RESULT =====

    % Report changes in .mat structure
    if isRaw
        DataMat.F = sFile;
    else
        DataMat.Events = sFile.events;
    end
    DataMat = rmfield(DataMat, 'Time');
    % Save file definition
    bst_save(file_fullpath(sInputs(iInput).FileName), DataMat, 'v6', 1);

end

OutputFiles = {sInputs.FileName};
end


function artifact_times = Compute(F, TimeVector, signal_type)
    if strcmp(signal_type, 'ECG')
        artifact_times = ecg_artifact_detection(F, TimeVector);
    elseif strcmp(signal_type, 'BP')
        artifact_times = bp_staircase_detection(F, TimeVector);
    end
end

function artifact_time = ecg_artifact_detection(ecg, time)
%ECG_ARTIFACT_DETECTION Detect artifacts (abrupt transients / saturation) in an ECG signal.
%
%   ARTIFACTS = ECG_ARTIFACT_DETECTION(ECG, T) analyzes the ECG signal ECG,
%   sampled at time points T (seconds), and returns a 2xM matrix ARTIFACTS
%   where each column is one detected artifact segment:
%       ARTIFACTS(1,i) = start time of artifact i
%       ARTIFACTS(2,i) = end time of artifact i
%
%   Detection combines two complementary criteria:
%     1) Time-domain clipping/saturation detection: flat runs of samples
%        sitting at the signal's extreme (rail) values.
%     2) Time-frequency broadband-energy detection using a continuous
%        wavelet transform with an analytic Morlet wavelet ('amor'):
%        segments where total scalogram energy is an outlier AND that
%        energy is spread outside the physiological QRS band (~8-25 Hz)
%        are flagged, since genuine QRS complexes concentrate energy in
%        that band while abrupt/saturation transients are broadband.
%
%   Requires the Wavelet Toolbox (uses CWT) and Signal Processing basics.
%
%   Example:
%       artifacts = ecg_artifact_detection(ecg, t);
%       for k = 1:size(artifacts,2)
%           fprintf('Artifact %d: %.3f s - %.3f s\n', k, artifacts(1,k), artifacts(2,k));
%       end

    % ---- Input checks ---------------------------------------------------
    ecg = ecg(:);
    time   = time(:);
    if numel(ecg) ~= numel(time)
        error('ecg_artifact_detection:sizeMismatch', ...
            'ecg and t must have the same number of samples.');
    end
    N = numel(ecg);
    if N < 10
        artifact_time = zeros(2,0);
        return;
    end

    fs = 1 / median(diff(time));

    % ---- Tunable parameters for ECG ---------------------------------------------
    qrsBandHz   = [8 25];   % physiological QRS energy band
    ratioThresh = 0.4;      % below this, energy judged "broadband"
    madK        = 20;        % MAD multiplier for energy-outlier threshold
    clipRunMin  = 100;        % min consecutive samples to call it clipping
    clipTol     = 1e-6;     % relative tolerance for "sitting at the rail"
    mergeGapSec = 0.05;     % merge artifacts separated by less than this
    minDurSec   = 0.1;     % discard artifacts shorter than this

    % ---- 1) Clipping / saturation detection ------------------------------
    sigMax = max(ecg);
    sigMin = min(ecg);
    rangeSig = sigMax - sigMin;
    if rangeSig == 0
        rangeSig = 1; % degenerate constant-signal edge case
    end
    atRail = (abs(ecg - sigMax) <= clipTol*rangeSig) | ...
             (abs(ecg - sigMin) <= clipTol*rangeSig);

    clipMask = false(N,1);
    i = 1;
    while i <= N
        if atRail(i)
            j = i;
            while j <= N && atRail(j)
                j = j + 1;
            end
            if (j - i) >= clipRunMin
                clipMask(i:j-1) = true;
            end
            i = j;
        else
            i = i + 1;
        end
    end

    % ---- 2) Time-frequency broadband-energy detection --------------------
    try
        [wt, f] = cwt(ecg, 'amor', fs);
    catch ME
        error('ecg_artifact_detection:cwtFailed', ...
            ['CWT computation failed (requires Wavelet Toolbox). ' ...
             'Original error: %s'], ME.message);
    end

    power = abs(wt).^2;                       % scalogram (frequency x time)
    totalPower = sum(power, 1) + eps;

    qrsIdx = f >= qrsBandHz(1) & f <= qrsBandHz(2);
    if ~any(qrsIdx)
        qrsIdx = true(size(f)); % fallback if band falls outside available scales
    end
    qrsPower = sum(power(qrsIdx, :), 1);
    ratio = qrsPower ./ totalPower;

    energyMedian = median(totalPower);
    energyMAD    = median(abs(totalPower - energyMedian)) + eps;
    isEnergyOutlier = totalPower > (energyMedian + madK*energyMAD);
    isBroadband     = ratio < ratioThresh;

    waveletMask = (isEnergyOutlier & isBroadband).'; % column vector

    % ---- Combine masks -----------------------------------------------------
    artifactMask = clipMask | waveletMask;

    % ---- Mask -> intervals, merge, filter -----------------------------------
    artifact_time = mask_to_intervals(artifactMask, time);
    artifact_time = merge_intervals(artifact_time, mergeGapSec);
    artifact_time = drop_short_intervals(artifact_time, minDurSec);

end

function intervals = mask_to_intervals(mask, t)
    d = diff([0; mask(:); 0]);
    startIdx = find(d == 1);
    endIdx   = find(d == -1) - 1;
    if isempty(startIdx)
        intervals = zeros(2,0);
        return;
    end
    intervals = [t(startIdx).'; t(endIdx).'];
end

function out = merge_intervals(in, gap)
    if isempty(in)
        out = in;
        return;
    end
    [starts, idx] = sort(in(1,:));
    ends = in(2, idx);

    outStarts = starts(1);
    outEnds   = ends(1);
    k = 1;
    for i = 2:numel(starts)
        if starts(i) - outEnds(k) <= gap
            outEnds(k) = max(outEnds(k), ends(i));
        else
            k = k + 1;
            outStarts(k) = starts(i);
            outEnds(k) = ends(i);
        end
    end
    out = [outStarts; outEnds];
end

function artifacts = bp_staircase_detection(bp, t)
%BP_STAIRCASE_DETECTION Detect staircase-like (frozen/quantized) artifacts in a
%blood pressure (BP) waveform.
%
%   ARTIFACTS = BP_STAIRCASE_DETECTION(BP, T) analyzes the BP signal BP,
%   sampled at time points T (seconds), and returns a 2xM matrix ARTIFACTS
%   where each column is one detected artifact segment:
%       ARTIFACTS(1,i) = start time of artifact i
%       ARTIFACTS(2,i) = end time of artifact i
%
%   Unlike a saturation/burst artifact, a staircase artifact is a LOSS of
%   pulsatility: the signal freezes onto discrete plateaus (e.g. a stuck
%   sample-and-hold, a dropped-and-padded communication link, or coarse
%   quantization), typically lasting a few seconds. Detection therefore
%   looks for the opposite signature of a burst detector:
%
%     1) Local flatness: sustained runs where sample-to-sample change is
%        near zero (the literal "steps" of the staircase), merged across
%        the brief small jumps that separate plateaus.
%     2) Pulse-band energy underflow: a continuous wavelet transform
%        (analytic Morlet) is used to track energy in the physiological
%        BP pulse band (~0.5-12 Hz); sustained energy well BELOW the
%        local baseline indicates the waveform has stopped pulsing.
%
%   Requires the Wavelet Toolbox (uses CWT).
%
%   Example:
%       artifacts = bp_staircase_detection(bp, t);
%       for k = 1:size(artifacts,2)
%           fprintf('Staircase %d: %.3f s - %.3f s\n', k, artifacts(1,k), artifacts(2,k));
%       end

    % ---- Input checks ---------------------------------------------------
    bp = bp(:);
    t  = t(:);
    if numel(bp) ~= numel(t)
        error('bp_staircase_detection:sizeMismatch', ...
            'bp and t must have the same number of samples.');
    end
    N = numel(bp);
    if N < 10
        artifacts = zeros(2,0);
        return;
    end

    fs = 1 / median(diff(t));

    % ---- Tunable parameters (staircase artifacts, ~3-10 s typical) --------
    pulseBandHz = [0.5 12];          % physiological BP pulse band
    madK        = 3;                 % MAD multiplier for energy UNDERflow
    flatWinSec  = 0.2;                % min run length to call a local plateau
    flatTol     = 0.1;               % relative tolerance for "samples ~equal"
    mergeGapSec = 0.7;                % bridge brief inter-plateau jumps
    minDurSec   = 0.25;                % below typical 3-10 s range -> filters normal brief plateaus
    extent_sec = 0.1;

    mov_avg_win_sec = 0.5;
    mov_avg_win_size =  mov_avg_win_sec * fs;
    bp = movmean(bp, mov_avg_win_size);

    % ---- 1) Local flatness detection (time domain) ------------------------
    % sigScale = median(abs(bp - median(bp))) + eps; % robust amplitude scale (MAD)

    d1 = abs(diff(bp));
    [flatMask_d1, flatSample_d1] = detect_flat(d1, flatTol, round(flatWinSec*fs));
    flatMask_d1 = [false; flatMask_d1];
    flatSample_d1 = [false; flatSample_d1];
    
    d1_smooth = movmean(d1, mov_avg_win_size/10);
    d2 = abs(diff(d1_smooth));
    [flatMask_d2, flatSample_d2] = detect_flat(d2, flatTol, round(flatWinSec*fs));
    flatMask_d2 = [false; false; flatMask_d2];
    flatSample_d2 = [false; false; flatSample_d2];

    flatMask = flatMask_d1 | flatMask_d2;
    flatSample = flatSample_d1 | flatSample_d2;

    % % ---- 2) Pulse-band energy underflow (time-frequency) -------------------
    % try
    %     [wt, f] = cwt(bp, 'amor', fs);
    % catch ME
    %     error('bp_staircase_detection:cwtFailed', ...
    %         ['CWT computation failed (requires Wavelet Toolbox). ' ...
    %          'Original error: %s'], ME.message);
    % end
    % 
    % power = abs(wt).^2;
    % pulseIdx = f >= pulseBandHz(1) & f <= pulseBandHz(2);
    % if ~any(pulseIdx)
    %     pulseIdx = true(size(f)); % fallback if band falls outside available scales
    % end
    % pulsePower = sum(power(pulseIdx, :), 1);
    % 
    % energyMedian = median(pulsePower);
    % energyMAD    = median(abs(pulsePower - energyMedian)) + eps;
    % lowThresh    = max(energyMedian - madK*energyMAD, 0);
    % isLowEnergy  = pulsePower < lowThresh;
    % 
    % lowEnergyMask = isLowEnergy.'; % column vector
    % 
    % % ---- Combine masks -------------------------------------------------------
    % % Either signature (sustained local flatness, or sustained pulse-band
    % energy underflow) is sufficient evidence of a frozen/staircase segment.

    artifactMask = flatMask; % | lowEnergyMask;

    % ---- Mask -> intervals, merge, filter -------------------------------------
    artifacts = mask_to_intervals(artifactMask, t);
    artifacts = merge_intervals(artifacts, mergeGapSec);
    artifacts = drop_short_intervals(artifacts, minDurSec);

    % Account for boundary trimming because of smoothing
    % + add some extent
    extent = mov_avg_win_sec/2;
    artifacts(1, :) = max(0, artifacts(1,:) - extent);
    artifacts(2, :) = min(max(t), artifacts(2,:) + extent);
    % hold on; dd1 = [0;d1]; plot(t, normalize(bp-mean(bp)), 'g'); plot(t, dd1, 'b'); plot(t(flatMask), dd1(flatMask), 'xr');
    % hold on; dd1 = [0;d1]; plot(t, normalize(bp-mean(bp)), 'g'); plot(t, dd1, 'b'); plot(t(flatSample_d1), dd1(flatSample_d1), 'xr');
    % hold on; dd2 = [0; 0; d2]; plot(t, normalize(bp-mean(bp)), 'g'); plot(t, dd2, 'b'); plot(t(flatSample_d2), dd2(flatSample_d2), 'xr');
end

function [flatMask, flatSample] = detect_flat(d, flatTol, flat_win_size)
    N = length(d);
    sigScale = prctile(d, 95); 
    flatSample = d <= flatTol * sigScale;   

    flatRunMin = max(1, flat_win_size);
    flatMask = false(N,1);
    i = 1;
    while i <= N
        if flatSample(i)
            j = i;
            while j <= N && flatSample(j)
                j = j + 1;
            end
            if (j - i) >= flatRunMin
                flatMask(i:j-1) = true;
            end
            i = j;
        else
            i = i + 1;
        end
    end
end

function out = drop_short_intervals(in, minDur)
    if isempty(in)
        out = in;
        return;
    end
    dur = in(2,:) - in(1,:);
    out = in(:, dur >= minDur);
end
