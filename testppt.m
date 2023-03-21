

% inject_hough_rmode_wn
tic

fmin=650;
fmax=900;
% fmaxmax=850; ci aggiorniamo 
dt=1/(2*(fmax-fmin));


Tobs=1024;
Tfft=4;
tt=0:dt:Tobs;
len=length(tt);

% inject_hough_rmode_wn
% % nomefile='L-L1_GWOSC_O3a_4KHZ_R1-1253933056-4096.hdf5';
%  nomefile='H-H1_LOSC_4_V2-1126257414-4096.hdf5';
% noise=h5read(nomefile,'/strain/Strain');
% N=length(noise);
% dt=h5readatt(nomefile,'/strain/Strain','Xspacing');
% t0=h5readatt(nomefile,'/strain/Strain','Xstart');
% tt=(0:(N-1))*dt;
% % 

fdotmax=-1/Tfft^2;
fdotmin=-1/(2*Tfft)^2;

h0=3E-23   
fdot0=fdotmin+rand*(fdotmax-fdotmin);
f0=700+floor(rand*150); %ho cambiato la frequenza (650:790)
n=7;

kn=abs(fdot0)/(f0)^n; %parametri 
x0=f0^(1-n);

[strain,ff]=generate_long_transient(f0,kn,h0,n,tt);%ho aggiunto (j) per come avevo definito l'array h0(j)
ffinj=ff-fmin;
ph1=mod(cumsum(ffinj).*[tt(2)-tt(1) diff(tt)],1)*2*pi;
strain=strain.*sin(ph1);

strainmax=max(abs(strain));
ampnois=(1.3*10^(-23))/sqrt(dt);%*max(strain)
noise=randn(1,len)*ampnois;
data=strain+noise;

lfft=Tfft/dt;
nover=lfft/2;

tfps_clean=tfr_ps(strain,dt,lfft,nover,1,3);%devo aggiungere un ,3 alla fine 
tfps_clean.iniy=fmin;%+1/Tfft;
mask0=tfps_clean.Z/mean(tfps_clean.Z(:));

image_plot(tfps_clean);xlabel('t(s)'),ylabel('f(Hz)'),title('r-mode signal spectrogram')

tfps_noisy=tfr_ps(data,dt,lfft,nover,1,3);
tfps_noisy.Z=tfps_noisy.Z/mean(tfps_noisy.Z(:));
tfps_noisy.iniy=fmin;%+1/Tfft;
SNR_tfps=SNR_image_vs_mask(tfps_noisy.Z,mask0);

soglia=2.5;
[nx,ny]=size(tfps_noisy.Z); %nx e ny sono le dimensioni delle mie immagini

%metodo generico
Mi=-0.01;%fdotmin*tfps_clean.dx/tfps_clean.dy; %min
Mf=fdotmax*tfps_noisy.dx/tfps_noisy.dy;%max
tfps_gen=apply_trifilter(tfps_noisy,Mi,Mf,1);
SNR_gen=SNR_image_vs_mask(tfps_gen.Z,mask0);

% metodo custom
dff=diff(ff); 
Mi=max(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
Mf=min(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
tfps_cust=apply_trifilter(tfps_noisy,Mi,Mf,1);
SNR_gen=SNR_image_vs_mask(tfps_gen.Z,mask0);
SNR_cust=SNR_image_vs_mask(tfps_cust.Z,mask0);
% 
% image_plot(tfps_noisy);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of noisy r-mode signal:   SNR=%.3f',SNR_tfps))
% image_plot(tfps_gen);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('"Generalised" filtered spectrogram:   SNR=%.3f, cut=%.3f',SNR_gen,cut))
% image_plot(tfps_cust);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('"Custom" filtered spectrogram:   SNR=%.3f',SNR_cust))
% 
a1=tfps_gen.Z(:);
a1=sort(a1);
ind=round(length(a1)*0.92);
soglia_gen=a1(ind);
a2=tfps_cust.Z(:);
a2=sort(a2);
ind=round(length(a2)*0.92);
soglia_cust=a2(ind);
[peaks_noisy,inds,times]=LAB_create_peakmap_from_wnimg(tfps_noisy,soglia,1,1);
[peaks_gen,inds,times]=LAB_create_peakmap_from_wnimg(tfps_gen,soglia_gen,1,1);
[peaks_cust,inds,times]=LAB_create_peakmap_from_wnimg(tfps_cust,soglia_cust,1,1);

plot_triplets(peaks_noisy(1,:),peaks_noisy(2,:),peaks_noisy(3,:))
plot_triplets(peaks_gen(1,:),peaks_gen(2,:),peaks_gen(3,:))
plot_triplets(peaks_cust(1,:),peaks_cust(2,:),peaks_cust(3,:))

% [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_noisy,Tfft,[fmin fmax],n,kn);
% [hfdf_gen,gridk_gen,gridx_gen] = LAB_hough_and_peak_noantsour_nogd(peaks_gen,Tfft,[fmin fmax],n,kn);
% [hfdf_cust,gridk_cust,gridx_cust] = LAB_hough_and_peak_noantsour_nogd(peaks_cust,Tfft,[fmin fmax],n,kn);
% plot_triplets(hfdf(:),gridx(:),gridk(:))

% hfdf=Pixel(hfdf,x0,kn,gridx,gridk);
% hfdf_gen=Pixel(hfdf_gen,x0,kn,gridx_gen,gridk_gen);
% hfdf_cust=Pixel(hfdf_cust,x0,kn,gridx_cust,gridk_cust);
% 
% figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title('Hough tranform of spectrogram');
% figure,uimagesc(gridx_gen,gridk_gen,hfdf_gen'),colorbar,axis xy, title(compose('Hough tranform of "Generalised" filtered spectrogram'));
% figure,uimagesc(gridx_cust,gridk_cust,hfdf_cust'),colorbar,axis xy, title('Hough tranform of "Custom" filtered spectrogram');


% triangular_filter=generate_triangle(nx,ny,-ny/(nx*Mi),-ny/(nx*Mf),cut)+0.001;
% image_plot(triangular_filter);title(compose('Triangular filter with cut=%.3f',cut))
% pixel=Pixel(hfdf,x0,kn,gridx,gridk);
% massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
% media=mean_sub_Hough(pixel);
% dev_std=deviazione_std(pixel);
% cr=(massimo-media)/dev_std
% 
% pixel=Pixel(hfdf_gen,x0,kn,gridx_gen,gridk_gen);
% massimo=massimo_intorno(hfdf_gen,x0,kn,gridx_gen,gridk_gen);
% media=mean_sub_Hough(pixel);
% dev_std=deviazione_std(pixel);
% crtrif=(massimo-media)/dev_std
% 
% pixel=Pixel(hfdf_cust,x0,kn,gridx_cust,gridk_cust);
% massimo=massimo_intorno(hfdf_cust,x0,kn,gridx_cust,gridk_cust);
% media=mean_sub_Hough(pixel);
% dev_std=deviazione_std(pixel);
% crcust=(massimo-media)/dev_std
% 
