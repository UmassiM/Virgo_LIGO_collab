 
nomefile='H-H1_LOSC_4_V2-1126257414-4096.hdf5';
data=h5read(nomefile,'/strain/Strain');
N=length(data);
dt=h5readatt(nomefile,'/strain/Strain','Xspacing');
t0=h5readatt(nomefile,'/strain/Strain','Xstart');
t=(0:N-1)*dt;

% Different resolution spectra
[data_pow,f_data]=powsp(data,dt,1,0,1);
figure;loglog(f_data,data_pow.^0.5),grid on,xlim([5 2000]),ylim([1e-25 1e-18])
xlabel('freq (Hz)'),ylabel('Strain (Hz^{-1/2})')

[data_pow10,f10_data]=powsp(data,dt,1,0,10);
figure;loglog(f10_data,data_pow10.^0.5),grid on,xlim([5 2000]),ylim([1e-25 1e-18])
xlabel('freq (Hz)'),ylabel('Strain (Hz^{-1/2})')

% Spettrogram creation
Tfft=4;
lfft=ceil(Tfft/dt);
lfft=floor(lfft/2)*2;
nover=lfft/2;
tfps_0=tfr_ps(data,dt,lfft,nover,1,5);
image_plot(tfps_0);xlabel('s'),ylabel('Hz'),title('Data spectrogram')

[pow_med,f_med]=powsp_median(data_pow,f_data,128,2);
figure;loglog(f_med,pow_med.^0.5),grid on,xlim([5 2000]),ylim([1e-25 1e-18])
xlabel('freq (Hz)'),ylabel('Strain (Hz^{-1/2})'),title('Median-evalued spectrum')


figure;loglog(f_data,data_pow.^0.5,'DisplayName','FFT Square'),grid
xlim([5 2000]),ylim([1e-25 1e-18]),xlabel('freq (Hz)'),ylabel('Strain (Hz^{-1/2})')
hold on,loglog(f_med,pow_med.^0.5,'DisplayName','Median Spectrum'),legend,hold off

lengratio=length(data_pow)/length(pow_med);
iii=floor((1:length(data_pow))/lengratio)+1;
iii(iii>length(pow_med))=length(pow_med);
pow_med=pow_med(:);
ratio=data_pow./pow_med(iii);
figure,semilogy(f_data,ratio),grid