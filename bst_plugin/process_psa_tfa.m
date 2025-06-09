function varargout = process_psa_tfa( varargin )

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
sProcess.Comment     = 'Transfer Function Analysis';
sProcess.FileTag     = 'TFA';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Frequency';
sProcess.Index       = 4; %0: not shown, >0: defines place in the list of processes
sProcess.Description = '';
sProcess.isSeparator = 0; % add a horizontal bar after the process in the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'matrix', 'matrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options
sProcess.options.ref_channel.Comment = 'Reference Channel: ';
sProcess.options.ref_channel.Type    = 'text';
sProcess.options.ref_channel.Value   = '';

sProcess.options.transfer_channels.Comment = 'Channels with transfer signal (comma-separated types or names): ';
sProcess.options.transfer_channels.Type    = 'text';
sProcess.options.transfer_channels.Value   = '';
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
    
    dt = diff(sDataIn.Time(1:2));
    channels = in_bst_channel(sInputs(iInput).ChannelFile);
    nb_channels = size(channels.Channel, 2);
        
    ref_sig_idx_chan = channel_find(channels.Channel, sProcess.options.ref_channel.Value);
    transfer_sig_idx_chan = channel_find(channels.Channel, sProcess.options.transfer_channels.Value);
    
    index_labels = tfa_index_labels();
    ref_signal = sDataIn.F(ref_sig_idx_chan, :);
    for ichan=1:length(transfer_sig_idx_chan)
        idx_chan = transfer_sig_idx_chan(ichan);
        transfer_signal = sDataIn.F(idx_chan, :);
        tfa_out = Compute(ref_signal, transfer_signal, 1/dt);
        chan_name = channels.Channel(idx_chan).Name;
        for ilabel=1:length(index_labels)
            index_label = index_labels{ilabel};
            tf_values.(protect_field_label([index_label '_' chan_name])) = tfa_out.(index_label);
        end
    end   
    
    % Save as matrix
    comment = [sInputs(iInput).Comment  ' | TFA' ];
    OutputFiles{iInput} = nst_save_table_in_bst(struct2table(tf_values), sInputs(iInput).SubjectName, sInputs(iInput).Condition,...
                                                  comment);
end

end

function labels = tfa_index_labels()
labels = {'Gain_vlf', 'Phase_vlf', 'Gain_lf', 'Phase_lf', 'Gain_hf', 'Phase_hf', ...
          'Gain_vlf_norm', 'Gain_lf_norm', 'Gain_hf_norm'};
end

function label = protect_field_label(label)
label = strrep(label, '.', '_dot_');
end



%% ===== Compute =====
function [tfa_out] = Compute(bp_signal, cbfi_signal, fs)

params.plot = 0;
tfa_out = tfa_car(bp_signal, cbfi_signal, fs, params);

end


function tfa_out = tfa_car(ABP,CBFV,fs,params)

% function tfa_out=psa_tfa_car(ABP,CBFV,fs,params)
% Transfer function analysis for cerebral autoregulation
% ABP and CBFV are in mmHg and cm/s, respecitvely, fs is the sampling
% frequency and params holds parameters for tfa analysis. Params is
% optional, and if not given will automatically use the default settings
% as given in the White Paper (2015).

% Prepared by David Simpson (ds@isvr.soton.ac.uk), Sept. 2015
% Adapted by Thomas Vincent, 2025

% default parameters taken from Consensus Paper
default_params.vlf=[0.02,0.07];
default_params.lf=[0.07,0.2];
default_params.hf=[0.2,0.5];
default_params.detrend=0;
default_params.spectral_smoothing=3; % triangular filter of this length
default_params.coherence2_thresholds=[3:15;0.51,0.40,0.34,0.29,0.25,0.22,0.20,0.18,0.17,0.15,0.14,0.13,0.12]';% from consensus paper
; % the first column refers to the number of 50% overlapping hanning windows, with 3-point spectral smoothing; only coherence values above this are used in calculating the mean gain and phase
default_params.apply_coherence2_threshold=1;
default_params.remove_negative_phase=1;
default_params.remove_negative_phase_f_cutoff=0.1;
default_params.normalize_ABP=0;
default_params.normalize_CBFV=0;
default_params.window_type='hanning';% alternatives 'Boxcar'
default_params.window_length=102.4;% in s
default_params.overlap=59.99;% overlap in % (when overlap_adjust is on, this is adjusted down). Use 59.99% rather than 60 so with with data corresponding to 5 windows with 50% overlap, 5 windows are chosen
default_params.overlap_adjust=1;
default_params.plot=0;
default_params.plot_f_range=[0,0.5];
default_params.plot_title='';

% coordinates for the plots in relative screen units
subplot1=[.15,0.65,0.75,0.22];
subplot2=subplot1+[0,-0.25,0,0];
subplot3=subplot2+[0,-0.25,0,0];

% set all missing parameters in to the default values
if nargin<4
    params=default_params;
end
if ~isfield(params,'vlf')
    params.vlf=default_params.vlf;
end
if ~isfield(params,'lf')
    params.lf=default_params.lf;
end
if ~isfield(params,'hf');
    params.hf=default_params.hf;
end
if ~isfield(params,'detrend');
    params.detrend=default_params.detrend;
end
if ~isfield(params,'spectral_smoothing');
    params.spectral_smoothing=default_params.spectral_smoothing;
end
if ~isfield(params,'coherence2_thresholds');
    params.coherence2_thresholds=default_params.coherence2_thresholds;
end
if ~isfield(params,'apply_coherence2_threshold');
    params.apply_coherence2_threshold=default_params.apply_coherence2_threshold;
end
if ~isfield(params,'remove_negative_phase');
    params.remove_negative_phase=default_params.remove_negative_phase;
end
if ~isfield(params,'remove_negative_phase_f_cutoff');
    params.remove_negative_phase_f_cutoff=default_params.remove_negative_phase_f_cutoff;
end
if ~isfield(params,'normalize_ABP');
    params.normalize_ABP=default_params.normalize_ABP;
end
if ~isfield(params,'normalize_CBFV');
    params.normalize_CBFV=default_params.normalize_CBFV;
end
if ~isfield(params,'window_type');
    params.window_type=default_params.window_type;
end
if ~isfield(params,'window_length');
    params.window_length=default_params.window_length;
end
if ~isfield(params,'overlap');
    params.overlap=default_params.overlap;
end
if ~isfield(params,'overlap_adjust');
    params.overlap_adjust=default_params.overlap_adjust;
end
if ~isfield(params,'default_params.plot');
    params.plot=default_params.plot;
end
if ~isfield(params,'plot_f_range');
    params.plot_f_range=default_params.plot_f_range;
end
if ~isfield(params,'plot_title');
    params.plot_title=default_params.plot_title;
end

tfa_out.Mean_abp=mean(ABP);
tfa_out.Std_abp=std(ABP);

if params.detrend
    ABP=detrend(ABP);
else
    ABP=ABP-mean(ABP);
end
if params.normalize_ABP==1
    ABP=(ABP/mean(ABP))*100;
end

tfa_out.Mean_cbfv=mean(CBFV);
tfa_out.Std_cbfv=std(CBFV);
if params.detrend
    CBFV=detrend(CBFV);
else
    CBFV=CBFV-mean(CBFV);
end
if params.normalize_CBFV==1
    CBFV=(CBFV/mean(CBFV))*100;
end

window_length=round(params.window_length*fs);
if strcmp(upper(params.window_type),'HANNING');
    wind=hanning_car(window_length);
end
if strcmp(upper(params.window_type),'BOXCAR');
    wind=boxcar(window_length);
end
overlap=params.overlap;
if params.overlap_adjust==1
    L=floor((length(ABP)-window_length)/(window_length*(1-params.overlap/100)))+1;
    if L>1
        shift=floor((length(ABP)-window_length)/(L-1));
        overlap=(window_length-shift)/window_length*100;
    end
end
tfa_out.overlap=overlap;

overlap=overlap/100;   
M_smooth=params.spectral_smoothing;
N_fft=window_length;
% keyboard


[H,C,f,Pxx,Pxy,Pyy,no_windows]=tfa1(ABP,CBFV,wind,overlap,M_smooth,fs,N_fft);

% keyboard



tfa_out.H=H;
tfa_out.C=C;
tfa_out.f=f;
tfa_out.Pxx=Pxx;
tfa_out.Pyy=Pyy;
tfa_out.Pxy=Pxy;
tfa_out.No_windows=no_windows;

% find the mean values in each frequency band
i=find(params.coherence2_thresholds(:,1)==no_windows);
if isempty(i)
   warning('No coherence threshold defined for the number of windows obtained - all frequencies will be included');
   coherence2_threshold=0;
else
   coherence2_threshold=params.coherence2_thresholds(i,2);
end

G=H; % save for plotting below
if params.apply_coherence2_threshold
i=find(abs(C).^2 < coherence2_threshold); % exclude low coherence
H(i)=nan;
% keyboard
end
P=angle(H);
if params.remove_negative_phase; % exclude negative phase below cut-off frequency
    n=find(f<params.remove_negative_phase_f_cutoff);
    k=find(P(n)<0);
    if ~isempty(k);
        P(n(k))=nan;
    end
end;
    
i=find(f>=params.vlf(1) & f<params.vlf(2));
tfa_out.Gain_vlf=mean(abs(H(i)), "omitnan");
tfa_out.Phase_vlf=mean(P(i), "omitnan")/(2*pi)*360;
tfa_out.Coh2_vlf=mean(abs(C(i)).^2, "omitnan");
tfa_out.P_abp_vlf=2*sum(Pxx(i))*f(2);
tfa_out.P_cbfv_vlf=2*sum(Pyy(i))*f(2);


i=find(f>=params.lf(1) & f<params.lf(2));
tfa_out.Gain_lf=mean(abs(H(i)), "omitnan");
tfa_out.Phase_lf=mean(P(i), "omitnan")/(2*pi)*360;
tfa_out.Coh2_lf=mean(abs(C(i)).^2, "omitnan");
tfa_out.P_abp_lf=2*sum(Pxx(i))*f(2);
tfa_out.P_cbfv_lf=2*sum(Pyy(i))*f(2);

i=find(f>=params.hf(1) & f<params.hf(2));
tfa_out.Gain_hf=mean(abs(H(i)), "omitnan");
dummy=angle(H(i));
tfa_out.Phase_hf=mean(P(i), "omitnan")/(2*pi)*360;
tfa_out.Coh2_hf=mean(abs(C(i)).^2, "omitnan");
tfa_out.P_abp_hf=2*sum(Pxx(i))*f(2);
tfa_out.P_cbfv_hf=2*sum(Pyy(i))*f(2);

if params.normalize_CBFV
    tfa_out.Gain_vlf_norm=tfa_out.Gain_vlf;
    tfa_out.Gain_lf_norm=tfa_out.Gain_lf;
    tfa_out.Gain_hf_norm=tfa_out.Gain_hf;
    tfa_out.Gain_vlf_not_norm=tfa_out.Gain_vlf*tfa_out.Mean_cbfv/100;
    tfa_out.Gain_lf_not_norm=tfa_out.Gain_lf*tfa_out.Mean_cbfv/100;
    tfa_out.Gain_hf_not_norm=tfa_out.Gain_hf*tfa_out.Mean_cbfv/100;
else
    tfa_out.Gain_vlf_not_norm=tfa_out.Gain_vlf;
    tfa_out.Gain_lf_not_norm=tfa_out.Gain_lf;
    tfa_out.Gain_hf_not_norm=tfa_out.Gain_hf;
    tfa_out.Gain_vlf_norm=tfa_out.Gain_vlf/tfa_out.Mean_cbfv*100;
    tfa_out.Gain_lf_norm=tfa_out.Gain_lf/tfa_out.Mean_cbfv*100;
    tfa_out.Gain_hf_norm=tfa_out.Gain_hf/tfa_out.Mean_cbfv*100;
end

    

if params.plot
    fig=figure;
    % set(fig,'Position',pos3);
    t=(0:length(ABP)-1)/fs;
    plot(t,ABP,t,CBFV,'r:');
    title(params.plot_title);
    xlabel('time (s)');
    legend('ABP','CBFV');
    axis tight
    fig=figure;
    ax(1)=subplot('position',subplot1);%(3,1,1);
    plot(f,abs(G));
    title(params.plot_title);
    ylabel('Gain');
    
    hold on
    plot(params.vlf,[1,1]*tfa_out.Gain_vlf,':r');
    plot(params.lf,[1,1]*tfa_out.Gain_lf,':r');
    plot(params.hf,[1,1]*tfa_out.Gain_hf,':r');
    set(ax(1),'XTickLabel','');
    axis tight
    ax(2)=subplot('position',subplot2);%(3,1,2);
    set(ax(2),'XTickLabel','');
    plot(f,angle(G)/(2*pi)*360);
    hold on
    ylabel('Phase (deg)');
    plot(params.vlf,[1,1]*tfa_out.Phase_vlf,':r');
    plot(params.lf,[1,1]*tfa_out.Phase_lf,':r');
    plot(params.hf,[1,1]*tfa_out.Phase_hf,':r');
    set(ax(2),'XTickLabel','');
    axis tight
    
    ax(3)=subplot('position',subplot3);%(3,1,3);
    plot(f,abs(C).^2);
    ylabel('|Coh|^2');
    hold on
    plot(params.vlf,[1,1]*tfa_out.Coh2_vlf,':r');
    plot(params.lf,[1,1]*tfa_out.Coh2_lf,':r');
    plot(params.hf,[1,1]*tfa_out.Coh2_hf,':r');
    plot(params.plot_f_range,coherence2_threshold*ones(1,2),'--k');
    axis tight
    
    xlabel('frequency(Hz)');
    linkaxes(ax,'x');
    xlim(params.plot_f_range);
    
end

end


function [H,C,f,Pxx,Pxy,Pyy,no_windows]=tfa1(x,y,wind,overlap,M_smooth,fs,Nfft);
% Transfer function analysis using the Welch method
% 
% [H,C,f,Pxx,Pxy,Pyy,no_windows]=tfa(x,y,T,overlap,fs,Nfft_opt);
% where x is the input signal, y the output signal, 
% M the window-length (in seconds, using a boxcar/rectangular window); if may be vector, and then will
% used as the window (e.g. hanning(256)), 
% overlap is the overlap (as a fraction of 1), 
% M_smooth - length of the triangular smoothing function. Must be odd
% number. 
% Nfft=0.. make the fft as long as the window (default value), else it defines the length of the fft (> window length)
% 
% H is the complex frequency response (transfer function)
% C is the coherence
% f is the corresponding frequency vector. 
% Pxx, Pxy, Pyy are the auto and cross-spectral densities. 
% no_windows is the number of windows used. 
% (C) David Simpson, December 2014.

if length(wind)==1
    M=wind;
    wind=boxcar(wind);
else
    M=length(wind);
end
    

if nargin<7 
    Nfft=0;
end
if Nfft==0
    Nfft=M;
end

[C,f,no_windows]=welch1(x,y,wind,overlap,fs);
Pxx=C.Pxx;
Pyy=C.Pyy;
Pxy=C.Pxy;

% [Pxy,f]=cpsd(y,x,wind,round(M*overlap),Nfft,fs,'twosided');
% [Pxx,f]=cpsd(x,x,wind,round(M*overlap),Nfft,fs,'twosided');
% [Pyy,f]=cpsd(y,y,wind,round(M*overlap),Nfft,fs,'twosided');

if M_smooth>1;
h=ones(floor((M_smooth+1)/2),1);
h=h/sum(h);
Pxx1=Pxx;
Pxx1(1)=Pxx(2);
Pyy1=Pyy;
Pyy1(1)=Pyy(2);
% Pxy1=abs(Pxy);
% Pxy1abs(1)=abs(Pxy(2));
Pxy1=(Pxy);
Pxy1(1)=Pxy(2);

Pxx1=filtfilt(h,1,Pxx1);
Pyy1=filtfilt(h,1,Pyy1);
Pxy1=filtfilt(h,1,Pxy1);

Pxx1(1)=Pxx(1);
Pxx=Pxx1;
Pyy1(1)=Pyy(1);
Pyy=Pyy1;
% Pxy1=Pxy1abs.*Pxy./abs(Pxy);
% Pxy1(1)=Pxy(1);
% Pxy1=Pxy1abs.*Pxy./abs(Pxy);
Pxy1(1)=Pxy(1);

% keyboard
Pxy=Pxy1;

end

% keyboard

H=Pxy./Pxx;
C=Pxy./(abs(Pxx.*Pyy).^0.5);

% no_windows=0;
% w_end=M;
% while w_end<=length(x);
%     no_windows=no_windows+1;
%     w_end=w_end+round(M-M*overlap);
% end

end

function [C,f,L]=welch1(x,y,window,overlap,fs,Nfft)
% function [C,f,L]=welch(x,y,window,overlap,fs);
% x and y are input signals
% window is either the window function (e.g. hanning(256)), or the length
% of the window - Boxcar is the default window in that case
% overlap is the fractional overlap, and fs the sampling frequency.
% PSC, CPSD and coherence are returned in C. L is the number of windows. 

if length(window)==1
    M=window;
    window=boxcar(M);
else
    M=length(window);
end
if nargin<6
    Nfft=M;
end

shift=round((1-overlap)*M);
x=x(:);
y=y(:);
window=window(:);
N=length(x);

X=fft(x(1:M).*window);
Y=fft(y(1:M).*window);
L=1;
if shift>0
i_start=1+shift;
while i_start+M-1 <= N
  X=[X,fft(x(i_start:i_start+M-1).*window,Nfft)];
  Y=[Y,fft(y(i_start:i_start+M-1).*window,Nfft)];
  i_start=i_start+shift;
  L=L+1;
end
end
f=(0:Nfft-1)'/Nfft*fs;
if L==1
C.Pxx=(X.*conj(X))/L/sum(window.^2)/fs;
C.Pyy=(Y.*conj(Y))/L/sum(window.^2)/fs;
C.Pxy=(conj(X).*Y)/L/sum(window.^2)/fs;
C.coh=C.Pxy./((abs((C.Pxx.*C.Pyy))).^0.5);
else
C.Pxx=sum(X.*conj(X),2)/L/sum(window.^2)/fs;
C.Pyy=sum(Y.*conj(Y),2)/L/sum(window.^2)/fs;
C.Pxy=sum(conj(X).*Y,2)/L/sum(window.^2)/fs;
C.coh=C.Pxy./((abs((C.Pxx.*C.Pyy))).^0.5);
end

end

function w=hanning_car(M)
% w=(1-cos(2*pi*([0:M-1]'+0.5)/M))/2;
w=(1-cos(2*pi*((0:M-1)')/M))/2;
end



function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end
