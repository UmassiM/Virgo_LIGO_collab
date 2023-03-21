function tfps=tfr_ps(y,dt,lfft,nover,res,win)
% TFR_PS  time-frequency power spectrum
%
%     tfr_ps(y,dt,lfft,nover,res,win)
%
%    y      input data
%    dt     sampling time
%    lfft   length of ffts in samples (even)
%    nover  number of samples for overlapping (0 -> lfft/2)
%    res    over-resolution (> 1; 1 = natural)
%    win    window  (0 -> no, 1 -> bartlett, 2 -> hanning (def), 3 -> flatcos, 4 -> tukey, 5 -> gauss)

% Sapienza Università di Roma
% Laboratorio di Segnali e Sistemi II
% Author: Sergio Frasca - 2018

if ~exist('win','var')
    win=5;
end

y=y(:);

lfft=floor(lfft/2+1)*2;

if nover == 0
    nover=lfft/2;
end

n=length(y);

Dt=nover*dt;

nt=2*round(0.5*Dt/dt);
nt2=nt/2;
lfft1=round(lfft*res);
nf=lfft1;
ir=0;
if isreal(y)
    nf=ceil(nf/2);
    ir=1;
end
Nt=floor(n/lfft)-1;
tfps0=zeros(Nt,nf);
dfr=1/(dt*lfft1);
w=ps_window(win,lfft);
i=0;
j=0;

while i <= n-lfft
    j=j+1;
    x=y(i+1:i+lfft);
    x=x.*w.';
    xm=mean(x)*ir;
    x=x-xm;
    x(lfft+1:lfft1)=0;
    x=fft(x);
    tfps0(j,:)=abs(x(1:nf)).^2*dt/nt;
    i=i+nover;
end

tfps.Z=tfps0;
tfps.dx=Dt;
tfps.dy=dfr;
tfps.inix=0;
tfps.iniy=0;