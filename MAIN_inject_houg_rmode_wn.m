tic

% inject_hough_rmode_wn
n=7;
Tfft=4;
fmin=650;
fmax=900
dt=1/(2*(fmax-fmin));
Tobs=1024;
Tfft=4;
% Time vector
tt=0:dt:Tobs;
len=length(tt);

h0=1E-23;

lfft=Tfft/dt;
nover=lfft/2;
fdotmax=-1/Tfft^2;
fdotmin=-1/(2*Tfft)^2;

fdot0=fdotmin+rand*(fdotmax-fdotmin);
f0=fmin+50+rand*150;
x0=f0^(1-n);
kn=abs(fdot0)/(f0)^n;
[strain,ff]=generate_long_transient(f0,kn,h0,n,tt);
ph1=mod(cumsum(ff).*[tt(2)-tt(1) diff(tt)],1)*2*pi;
strain=strain.*sin(ph1);

% Simulation of noise choosing its amplitude with respect to signal
strainmax=max(abs(strain));
ampnois=50*max(strain);
noise=randn(1,len)*ampnois;

noismean=mean(noise);
noisstd=std(noise);
% amp_order=floor(log10(ampnois));
% ampord=10^amp_order;
% figure,histogram(noise/ampord),legend(compose('Mean=%.3e   Std=%.3e',[noismean noisstd])),grid
% title('Histogram of generated noise'),xlabel(compose('noise strain (e%d)',amp_order)),ylabel('counts')
SNR_strain=(strainmax-noismean)/noisstd;

% Injection of signal into noise and plot
data=strain+noise;
% figure,plot(tt,data),grid,hold on,plot(tt,strain,'r'),hold off
% xlabel('t(s)'),title(compose('r-mode signal with noise:   SNR=%.3f',SNR_strain))
lfft=Tfft/dt;
nover=lfft/2;
% Computation of spectrogram of the signal alone
tfps_clean=tfr_ps(strain,dt,lfft,nover,1);
% image_plot(tfps_clean);xlabel('t(s)'),ylabel('f(Hz)'),title('r-mode signal spectrogram')
% Creation of a mask, needed to calculate the image SNR
mask0=tfps_clean.Z/mean(tfps_clean.Z(:));
% Computation of spectrogram of full data, SNR calculation, plot
tfps_noisy=tfr_ps(data,dt,lfft,nover,1);
tfps_noisy.Z=tfps_noisy.Z/mean(tfps_noisy.Z(:));
SNR_tfps_0=SNR_image_vs_mask(tfps_noisy.Z,mask0);
image_plot(tfps_noisy);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of r-mode signal with noise:   SNR=%.3f',SNR_tfps_0))

% tfps_r=tfr_ps_r(data,dt,lfft,nover,1);
% tfps_r.Z(tfps_r.Z>20)=20;
% SNR_tfps_0=SNR_image_vs_mask(tfps_r.Z,mask0);
% image_plot(tfps_r);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of r-mode signal with noise:   SNR=%.3f',SNR_tfps_0))

[nx,ny]=size(tfps_clean.Z);
Mi=-0.01;
Mf=fdotmax*tfps_noisy.dx/tfps_noisy.dy;
triangular_filter=generate_triangle(nx,ny,-ny/(nx*Mi),-ny/(nx*Mf),1)+0.001;
tfps_trif=apply_trifilter(tfps_noisy,Mi,Mf,1);
SNR_tfps_trif=SNR_image_vs_mask(tfps_trif.Z,mask0);
enhance_trif=SNR_tfps_trif/SNR_tfps_0;
image_plot(tfps_trif);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Triangular-filtered image:   SNR=%.3f   enhance=%.3f',[SNR_tfps_trif enhance_trif]))

%distribuzioni 
% aaa=tfps_noisy.Z(:);
% figure,histogram(aaa)
% bbb=tfps_trif.Z(:);
% figure,histogram(bbb)

% Only the part of the map within the frequency band is selected
% i_fmin=round((fmin-tfps_noisy.iniy)/tfps_noisy.dy)+1;
% i_fmax=round((fmax-tfps_noisy.iniy)/tfps_noisy.dy)+1;
% tfps_noisy.Z=tfps_noisy.Z(:,i_fmin:i_fmax);
tfps_noisy.iniy=fmin;
% Peakmap creation
[peaks,inds,times]=LAB_create_peakmap_from_wnimg(tfps_noisy,2.5,1,1);
[hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks,Tfft,[fmin fmax],n,kn);
figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title(compose('ratio'));

pixel=Pixel(hfdf,x0,kn,gridx,gridk);
massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
media=mean_sub_Hough(pixel);
dev_std=deviazione_std(pixel);
cr=(massimo-media)/dev_std

a=tfps_trif.Z(:);
a=sort(a);
ind=round(length(a)*0.92);
thresh=a(ind);
[peaks,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,thresh,1,1);
[hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks,Tfft,[fmin fmax],n,kn); 
figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title(compose('ratio'));

