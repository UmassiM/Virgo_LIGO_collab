function [ps,frs,nfft]=powsp(in,dt,res,win,npiec,pref)
% simple power spectrum
%
%   [ps,frs,nfft]=powsp(in,dt,res,win,npiec,pref)
%                                                     
%   in     inpt vector
%   dt     sampling time
%   res    resolution ? 1
%   win    window (0 no, 1 bartlett, 2 hann, 3 flatcos, 5 gauss)
%   npiec  number of pieces (def 1)
%   pref   prefiltering (0 no, 1 no mean (def), 2 no trend)

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

if ~exist('win','var')
    win=2;
end
if ~exist('npiec','var')
    npiec=1;
end
if ~exist('pref','var')
    pref=1;
end
switch pref
    case 1
        in=in-mean(in);
    case 2
        in=detrend(in);
end
icreal=0;
if isreal(in)
    icreal=1;
end
n=length(in);
nn=floor(n/npiec);
in=in(:);
bias=mean(in);

nfft=round(nn*res);
dfr=1/(nfft*dt);
frs=(0:nfft-1)*dfr;

w=ps_window(win,nn);

w=w';
wnorm=sum(w.^2);
i1=1;
ps=zeros(nfft,1);
        
for i = 1:npiec
    in1=in(i1:i1+nn-1).*w; 
    in1(nn+1:nfft)=bias;
    ps0=abs(fft(in1)).^2;%size(in1),size(w)
    ps=ps+ps0;
    i1=i1+nn;
end

if icreal
    ps=ps(1:floor((nfft+1)/2));
    frs=frs(1:floor((nfft+1)/2));
end

ps=ps*dt/(npiec*wnorm);