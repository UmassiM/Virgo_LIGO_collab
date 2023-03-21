tic 

% inject_hough_rmode_wn
nomefile='H-H1_LOSC_4_V2-1126257414-4096.hdf5';
noise=h5read(nomefile,'/strain/Strain');
N=length(noise);
var.dt=h5readatt(nomefile,'/strain/Strain','Xspacing');
t0=h5readatt(nomefile,'/strain/Strain','Xstart');
var.tt=(0:N-1)*var.dt;

%parametri liberi
h0=3E-23; %strain

f0=700;
var.Tfft=4;
var.fmin=650;
var.fmax=900;
var.lfft=var.Tfft/var.dt;
var.nover=var.lfft/2;
var.fdotmax=-1/var.Tfft^2;
var.fdotmin=-1/(2*var.Tfft)^2;
fdot0=var.fdotmin+rand*(var.fdotmax-var.fdotmin);
var.f0=f0+floor(rand*150);
n=7;
var.x0=var.f0^(1-n);
var.n=7;
var.kn=abs(fdot0)/(var.f0)^var.n;
[strain,var.ff]=generate_long_transient(var.f0,var.kn,h0,var.n,var.tt);
ph1=mod(cumsum(var.ff).*[var.tt(2)-var.tt(1) diff(var.tt)],1)*2*pi;
strain=strain.*sin(ph1);

data=strain+noise';

% tfps=tfr_ps(data,var.dt,var.lfft,var.nover,1,3);
% mage_plot(tfps);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of Dataset 1'))
tfps_clean=tfr_ps(strain,var.dt,var.lfft,var.nover,1,3);
tfps_r=tfr_ps_r(data,var.dt,var.lfft,var.nover,1,3);
var.mask0=tfps_clean.Z/mean(tfps_clean.Z(:));

% image_plot(tfps_r);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of Dataset 1'))
tfps_c=tfps_r;
tfps_c.Z(tfps_c.Z>20)=20;
image_plot(tfps_c);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of Dataset 1'))

% taglia(var, tfps_r, 20, 20)

toc