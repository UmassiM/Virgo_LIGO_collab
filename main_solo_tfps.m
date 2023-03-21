tfps_clean=tfr_ps(strain,dtnew,lfft,nover,1);

tfps_noisy=tfr_ps(data,dtnew,lfft,nover,1);
tfps_noisy.iniy=fmin;
tfps_noisy.Z=tfps_noisy.Z/mean(tfps_noisy.Z(:));
mask0=tfps_clean.Z/mean(tfps_clean.Z(:));
SNR_tfps_0=SNR_image_vs_mask(tfps_noisy.Z,mask0);

tfps_noisy_r=tfr_ps_r(data,dtnew,lfft,nover,1);
tfps_noisy_r.iniy=fmin;
SNR_tfps_R=SNR_image_vs_mask(tfps_noisy_r.Z,mask0);

[peaks_noisy,inds,times]=LAB_create_peakmap_from_wnimg(tfps_noisy,2.5,1,1);
[hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_noisy,Tfft,[fmin fmax],n,kn);

[peaks_noisy_r,inds_r,times]=LAB_create_peakmap_from_wnimg(tfps_noisy_r,2.5,1,1);
[hfdf_r,gridk_r,gridx_r] = LAB_hough_and_peak_noantsour_nogd(peaks_noisy_r,Tfft,[fmin fmax],n,kn);

image_plot(tfps_noisy);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of r-mode signal:   SNR=%.3f',SNR_tfps_0))
image_plot(tfps_noisy_r);xlabel('t(s)'),ylabel('f(Hz)'),title(compose('Spectrogram of r-mode CLEANED signal:   SNR=%.3f',SNR_tfps_R))
figure,uimagesc(gridx,gridk,hfdf'),colorbar,axis xy, title(compose('normale'))
figure,uimagesc(gridx_r,gridk_r,hfdf_r'),colorbar,axis xy, title(compose('ratio'));
