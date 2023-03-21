tic

nomefile='H-H1_LOSC_4_V2-1126257414-4096.hdf5';
noise=h5read(nomefile,'/strain/Strain');
N=length(noise);
dt=h5readatt(nomefile,'/strain/Strain','Xspacing');
t0=h5readatt(nomefile,'/strain/Strain','Xstart');
tt=(0:N-1)*dt;

%parametri del segnale in esame
n=7;
Tfft=4;
fmin=650;
fmax=900;
lfft=Tfft/dt;
nover=lfft/2;
fdotmax=-1/Tfft^2;
fdotmin=-1/(2*Tfft)^2;

%parametri liberi
par.h_i=5E-25; %strain iniziale
par.h_f=7E-23;  %strain finale
h_r=par.h_f/par.h_i; %serve per l'incremento

par.nums=15; %quanti valori da testare in ind (riga 84)
par.iter=30; %numero di h0 da analizzare
par.runs=10; %numero di ripetizioni per ogni h0
par.ndiv=8; %in quanti clusters si vuole dividere il rumore. ||METTERE UNA POTENZA DI 2||
par.total=par.num*(par.ndiv-1)*par.runs*par.iter; %numero totale di tierazioni
par.cut=20; %soglia di taglio
par.set=20; %normalizzazione

t = 1:1:par.iter;
x = @(u) 10.^(tanh((u-1.5*mean(u))*1.5/par.iter)); 
y = @(u) par.h_f*u;
h0=par.h_i+y(x(t))-min(y(x(t))); %così scelgo strains più densi agli estremi
% semilogx(h0, crvech0,'.','MarkerSize',17)
% grid on

l=length(tt);
m=max(tt);
tobs=1024; 
r=floor(l/m*tobs); %divide tt in 4 
step=l/par.ndiv;
tt=tt(1:r); %per arrivare a t=tobs
crmat_trif=zeros(15,par.iter,par.runs,par.ndiv-1); %in questa quelli filtrati
crvec=zeros(nums,par.iter);
crstd=zeros(nums,par.iter);

for j=1:1:par.iter
    for i=(1:1:par.ndiv-1)
        subnoise=noise(1+(i-1)*step:r+(i-1)*step);
        for k=(1:1:par.runs)
            fdot0=fdotmin+rand*(fdotmax-fdotmin);
            f0=fmin+50+rand*150;
            x0=f0^(1-n);
            kn=abs(fdot0)/(f0)^n;
            [strain,ff]=generate_long_transient(f0,kn,h0(j),n,tt);
            ph1=mod(cumsum(ff).*[tt(2)-tt(1) diff(tt)],1)*2*pi;
            strain=strain.*sin(ph1);
            data=strain+subnoise';
            
            tfps_noisy=tfr_ps_r(data,dt,lfft,nover,1,3);
            i_fmin=round((fmin-tfps_noisy.iniy)/tfps_noisy.dy)+1;
            i_fmax=round((fmax-tfps_noisy.iniy)/tfps_noisy.dy)+1;
            tfps_noisy.Z=tfps_noisy.Z(:,i_fmin:i_fmax);
            tfps_noisy.Z(tfps_noisy.Z>par.cut)=par.set;
            tfps_noisy.iniy=fmin;
             
            [nx,ny]=size(tfps_noisy.Z);
            
%             dff=diff(ff);
%             Mi=max(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
%             Mf=min(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
%             
            Mi=-0.01; %piccolo valore diverso da 0 come derivata minima
            Mf=fdotmax*tfps_noisy.dx/tfps_noisy.dy;

            tfps_trif=apply_trifilter(tfps_noisy,Mi,Mf,1);
            a=tfps_trif.Z(:);
            a=sort(a);
            
            for p=0:1:par.num-1
                ind=round(length(a)*0.85+p/100);
                thresh=a(ind);
                [peaks,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,thresh,1,1);
                [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks,Tfft,[fmin fmax],n,kn);
                pixel=Pixel(hfdf,x0,kn,gridx,gridk);
                massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
                media=mean_sub_Hough(pixel);
                dev_std=deviazione_std(pixel);

                crmat_trif(p+1,j,k,i)=(massimo-media)/dev_std;
            end
        end 
    end
end

for p=1:1:par.num
    a=reshape(crmat_trif(1,:,:,:),[par.iter,par.runs,par.ndiv-1]);
    crvec(p,:)=mean(a,[2 3]);
    crstd(p,:)=std(a,[],[2 3]);
end 


Legend=cell(par.num,1);
figure,semilogx(h0, crvec(p,:),'.-','MarkerSize',10)
hold on
for p=2:1:15
    semilogx(h0, crvec(p,:),'.-','MarkerSize',10)
    Legend{p}=strcat('thresh=', num2str(0.85+(p-1)/100));
end
legend(Legend)
grid on

par.time=toc

%salvo le variabili utili (potrei anche salvare crmat) e le nomina con la
%data, per evitare sovrascritture con lo stesso nome
nome=datestr(now,'mm-dd HH-MM'); 
save(nome, 'h0','crmat','crmat_trif','par');
