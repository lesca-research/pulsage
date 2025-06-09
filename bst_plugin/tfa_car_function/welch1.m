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
f=[0:Nfft-1]'/Nfft*fs;
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



  

