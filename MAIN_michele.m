tic 

cut=20;
set=20;

% inject_hough_rmode_wn
nomefile='H-H1_LOSC_4_V2-1126257414-4096.hdf5';
noise=h5read(nomefile,'/strain/Strain');
N=length(noise);
dt=h5readatt(nomefile,'/strain/Strain','Xspacing');
t0=h5readatt(nomefile,'/strain/Strain','Xstart');
tt=(0:(N-1))*dt;
% 
a=floor(N/4);
noise=noise(1:a);
tt=tt(1:a);
%parametri liberi

h0=3E-23; %strain

n=7;
Tfft=4;
fmin=650;
fmax=900;
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

data=strain+noise';

tfps_noisy=tfr_ps_r(strain,dt,lfft,nover,1,3);
tfps_noisy.Z(tfps_noisy.Z>cut)=set;
image_plot(tfps_noisy);xlabel('t(s)'),ylabel('f(Hz)')

% i_fmin=round((fmin-tfps_noisy.iniy)/tfps_noisy.dy)+1;
% i_fmax=round((fmax-tfps_noisy.iniy)/tfps_noisy.dy)+1;
% tfps_noisy.Z=tfps_noisy.Z(:,i_fmin:i_fmax);
% tfps_noisy.iniy=fmin;
% 
% [peaks_noisy,inds,times]=LAB_create_peakmap_from_wnimg(tfps_noisy,2.5,1,1);
% [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_noisy,Tfft,[fmin fmax],n,kn);
% 
% figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title(compose('ratio'));
% 
% intorno=Intorno(x0,kn,gridx,gridk,hfdf);
% pixel=Pixel(hfdf,x0,kn,gridx,gridk);
% massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
% media=mean_sub_Hough(pixel);
% mediana=mediana_sub_Hough(pixel);
% dev_std=deviazione_std(pixel);
% cr=(massimo-media)/dev_std
% 
% [nx,ny]=size(tfps_noisy.Z);
% 
% % Mi=-0.01; %piccolo valore diverso da 0 come derivata minima
% % Mf=fdotmax*tfps_noisy.dx/tfps_noisy.dy;
% dff=diff(ff);
% Mi=max(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
% Mf=min(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
%             
% triangular_filter=generate_triangle(nx,ny,-ny/(nx*Mi),-ny/(nx*Mf),1)+0.001;
% tfps_trif=apply_trifilter(tfps_noisy,Mi,Mf,1);
% 
% a=tfps_trif.Z(:);
% a=sort(a);
% ind=round(length(a)*0.92);
% thresh=a(ind);
% [peaks_trif,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,thresh,1,1);
% [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_trif,Tfft,[fmin fmax],n,kn);
% %             
% 
% figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title(compose('ratio'));
% 
% intorno=Intorno(x0,kn,gridx,gridk,hfdf);
% pixel=Pixel(hfdf,x0,kn,gridx,gridk);
% massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
% media=mean_sub_Hough(pixel);
% mediana=mediana_sub_Hough(pixel);
% dev_std=deviazione_std(pixel);
% crtri=(massimo-media)/dev_std

toc