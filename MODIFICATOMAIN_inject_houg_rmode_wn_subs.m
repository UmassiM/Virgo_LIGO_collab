% inject_hough_rmode_wn
fsamp=4096;
dt=1/(2*fsamp);

fmin=600;
fmax=800;
dtnew=1/(2*(fmax-fmin));

Tobs=1000;
Tfft=4;
tt=0:dtnew:Tobs;
[a,len]=size(tt);

fdotmax=-1/Tfft^2;
fdotmin=-1/(2*Tfft)^2;
fdot0=fdotmin+rand*(fdotmax-fdotmin);
f0=750;
f0inj=f0-fmin;
DL=40;
n=7;
kn=abs(fdot0)/(f0)^n;
h0=1E-24;

[strain,ff]=generate_long_transient(f0,kn,h0,n,tt);
ffinj=ff-fmin;
ph1=mod(cumsum(ffinj).*[tt(2)-tt(1) diff(tt)],1)*2*pi;
strain=strain.*sin(ph1);

strainmax=max(abs(strain));
dsfact=fsamp/(fmax-fmin);
ampnois=50*max(strain);
noise=randn(1,len)*ampnois/sqrt(dsfact);
noismean=mean(noise);
noisstd=std(noise);
amp_order=floor(log10(ampnois));
ampord=10^amp_order;
SNR_strain=(strainmax-noismean)/noisstd;

data=strain+noise;

[strain_pow,f_strain]=powsp(strain,dtnew,1,0,1);
[noise_pow,f_noise]=powsp(noise,dtnew,1,0,1000);
f_strain=f_strain+fmin;%+1/Tfft;
f_noise=f_noise+fmin;%+1/Tfft;

lfft=Tfft/dtnew;
nover=lfft/2;

tfps_clean=tfr_ps(strain,dtnew,lfft,nover,1);
tfps_clean.iniy=fmin;%+1/Tfft;

mask0=tfps_clean.Z/mean(tfps_clean.Z(:));

tfps_noisy=tfr_ps(data,dtnew,lfft,nover,1);
tfps_noisy.Z=tfps_noisy.Z/mean(tfps_noisy.Z(:));
tfps_noisy.iniy=fmin;%+1/Tfft;
SNR_tfps_0=SNR_image_vs_mask(tfps_noisy.Z,mask0);
image_plot(tfps_noisy);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of r-mode signal with noise:   SNR=%.3f',SNR_tfps_0))


[peaks_noisy,inds,times]=LAB_create_peakmap_from_wnimg(tfps_noisy.Z,2.5,1,1);

[hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_noisy,Tfft,[fmin fmax],n,kn);
figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy
% image_gd2_nonunigrid(hfdf,gridk);

[nx,ny]=size(tfps_clean.Z);
Mi=0;%fdotmin*tfps_clean.dx/tfps_clean.dy;
Mf=fdotmax*tfps_clean.dx/tfps_clean.dy;
triangular_filter=generate_triangle(nx,ny,-ny/(nx*Mi),-ny/(nx*Mf),1)+0.001;
tfps_trif=apply_trifilter(tfps_noisy,Mi,Mf,1);
SNR_tfps_trif=SNR_image_vs_mask(tfps_trif.Z,mask0);
enhance_trif=SNR_tfps_trif/SNR_tfps_0;
image_plot(tfps_trif);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Triangular-filtered image:   SNR=%.3f   enhance=%.3f',[SNR_tfps_trif enhance_trif]))