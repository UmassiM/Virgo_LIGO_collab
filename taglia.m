function[]=taglia(vars, tfps_ratio, cut, set) 
%
%La funzione prende in input i parametri del segnale e lo spettrogramma
%(ratio) per farne hough e spettrogramma filtrato in due modi.
%
%cut è la soglia sopra la quale i valori di tfps.Z (ratio) vengono
%normalizzati a set
%
%fmin ed fmax servono per tagliare una finestra dello spettrogramma sulla
%quale queste operazioni. Per prendere tutto lo spettrogramma basta non
%inserire fmin ed fmax).
%
tic

i_fmin=round((vars.fmin-tfps_ratio.iniy)/tfps_ratio.dy)+1;
i_fmax=round((vars.fmax-tfps_ratio.iniy)/tfps_ratio.dy)+1;
vars.mask0=vars.mask0(:,i_fmin:i_fmax);
tfps_ratio.Z=tfps_ratio.Z(:,i_fmin:i_fmax);
tfps_ratio.iniy=vars.fmin;

tfps_ratio.Z(tfps_ratio.Z>cut)=set;

[peaks_ratio,inds,times]=LAB_create_peakmap_from_wnimg(tfps_ratio,2.5,1,1);
[hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_ratio,vars.Tfft,[vars.fmin vars.fmax],vars.n,vars.kn);
pixel=Pixel(hfdf,vars.x0,vars.kn,gridx,gridk);
massimo=massimo_intorno(hfdf,vars.x0,vars.kn,gridx,gridk);
media=mean_sub_Hough(pixel);
dev_std=deviazione_std(pixel);
cr=(massimo-media)/dev_std

SNR_tfps=SNR_image_vs_mask(tfps_ratio.Z,vars.mask0);
image_plot(tfps_ratio);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of Dataset 1: SNR=%.3f',SNR_tfps))
figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title(compose('Hough transform of unfiltered spectrogram'));

[nx,ny]=size(tfps_ratio.Z);

%filtro "tailored"
dff=diff(vars.ff);
Mi=min(dff)/vars.dt*tfps_ratio.dx/tfps_ratio.dy;
Mf=max(dff)/vars.dt*tfps_ratio.dx/tfps_ratio.dy;
triangular_filter=generate_triangle(nx,ny,-ny/(nx*Mi),-ny/(nx*Mf),1)+0.001;
tfps_trif=apply_trifilter(tfps_ratio,Mi,Mf,1);
a1=tfps_trif.Z(:);
a1=sort(a1);
ind=round(length(a1)*0.92);
soglia_trif=a1(ind);
[peaks_trif,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,soglia_trif,1,1);
[hfdf_trif,gridk_trif,gridx_trif] = LAB_hough_and_peak_noantsour_nogd(peaks_trif,vars.Tfft,[vars.fmin vars.fmax],vars.n,vars.kn);

pixel=Pixel(hfdf_trif,vars.x0,vars.kn,gridx_trif,gridk_trif);
massimo=massimo_intorno(hfdf_trif,vars.x0,vars.kn,gridx_trif,gridk_trif);
media=mean_sub_Hough(pixel);
dev_std=deviazione_std(pixel);
crcust=(massimo-media)/dev_std

SNR_tfps=SNR_image_vs_mask(tfps_trif.Z,vars.mask0);
image_plot(tfps_trif);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('"Custom" filtered spectrogram: SNR=%.3f',SNR_tfps))

%filtro generico
Mi=-0.01; %piccolo valore diverso da 0 come derivata minima
Mf=vars.fdotmax*tfps_ratio.dx/tfps_ratio.dy;
triangular_filter=generate_triangle(nx,ny,-ny/(nx*Mi),-ny/(nx*Mf),1)+0.001;
tfps_trif=apply_trifilter(tfps_ratio,Mi,Mf,1);
a1=tfps_trif.Z(:);
a1=sort(a1);
ind=round(length(a1)*0.92);
soglia_trif=a1(ind);
[peaks_trif,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,soglia_trif,1,1);
[hfdf_trif,gridk_trif,gridx_trif] = LAB_hough_and_peak_noantsour_nogd(peaks_trif,vars.Tfft,[vars.fmin vars.fmax],vars.n,vars.kn);

pixel=Pixel(hfdf_trif,vars.x0,vars.kn,gridx_trif,gridk_trif);
massimo=massimo_intorno(hfdf_trif,vars.x0,vars.kn,gridx_trif,gridk_trif);
media=mean_sub_Hough(pixel);
dev_std=deviazione_std(pixel);
crgen=(massimo-media)/dev_std

SNR_tfps=SNR_image_vs_mask(tfps_trif.Z,vars.mask0);
image_plot(tfps_trif);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('"Generalised" filtered spectrogram: SNR=%.3f',SNR_tfps))


toc
end
