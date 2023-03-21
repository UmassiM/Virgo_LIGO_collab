tic

%parametri liberi
par.h_i=5E-25; %strain iniziale
par.h_f=2E-23;  %strain finale
h_r=par.h_f/par.h_i; %serve per l'incremento

par.iter=10; %numero di h0 da analizzare
par.runs=10; %numero di ripetizioni per ogni h0
par.ndiv=8; %in quanti clusters si vuole dividere il rumore. ||METTERE UNA POTENZA DI 2||
par.total=(par.ndiv-1)*par.runs*par.iter; %numero totale di tierazioni
par.cut=20; %soglia di taglio
par.set=20; %normalizzazione
par.custom=0.92; %parametro per filtro adattato
par.generic=0.92; %parametro per filtro generico

%funzione che genera il vettore di h0 per ottenere una densità di punti efficace
t = 1:1:par.iter;
x = @(u) 10.^(tanh((u-1.5*mean(u))*1.5/par.iter)); 
y = @(u) par.h_f*u;
h0=par.h_i+y(x(t))-min(y(x(t))); %così scelgo strains più densi agli estremi
    
theFiles=dir('*.hdf5');
for p = 12
    par.nomefile = theFiles(p).name
    noise=h5read(par.nomefile,'/strain/Strain');
    N=length(noise);
    dt=h5readatt(par.nomefile,'/strain/Strain','Xspacing');
    t0=h5readatt(par.nomefile,'/strain/Strain','Xstart');
    tt=(0:N-1)*dt;
    l=length(tt);
    m=max(tt);
    tobs=1024; 
    r=floor(l/m*tobs); %divide tt in 4 
    step=l/par.ndiv;
    tt=tt(1:r); %per arrivare a t=tobs
    n=7;
    Tfft=4;
    fmin=650;
    fmax=900;
    lfft=Tfft/dt;
    nover=lfft/2;
    fdotmax=-1/Tfft^2;
    fdotmin=-1/(2*Tfft)^2;
    
    crmat=zeros(par.iter,par.runs,par.ndiv-1); %vengono salvati tutti i valori critici in questa matrice
    crmat_dff=zeros(par.iter,par.runs,par.ndiv-1); %in questa quelli filtrati
    crmat_gen=zeros(par.iter,par.runs,par.ndiv-1); %per un altro approccio al filtro
    
    for j=1:1:par.iter
        for i=(1:1:par.ndiv-1)
            %selezione del settore di rumore reale
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

                %spettrogramma
                tfps_noisy=tfr_ps_r(data,dt,lfft,nover,1,3);
                i_fmin=round((fmin-tfps_noisy.iniy)/tfps_noisy.dy)+1;
                i_fmax=round((fmax-tfps_noisy.iniy)/tfps_noisy.dy)+1;
                tfps_noisy.Z=tfps_noisy.Z(:,i_fmin:i_fmax);
                tfps_noisy.Z(tfps_noisy.Z>par.cut)=par.set;
                tfps_noisy.iniy=fmin;
                %hough e valore critico
                [peaks_noisy,inds,times]=LAB_create_peakmap_from_wnimg(tfps_noisy,2.5,1,1);
                [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks_noisy,Tfft,[fmin fmax],n,kn);
                pixel=Pixel(hfdf,x0,kn,gridx,gridk);
                massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
                media=mean_sub_Hough(pixel);
                dev_std=deviazione_std(pixel);
                crmat(j,k,i)=(massimo-media)/dev_std;

                %filtri triangolari
                [nx,ny]=size(tfps_noisy.Z);
                %filtro adattato (derivata discreta delle frequenze da strain generato
                dff=diff(ff);
                Mi=max(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
                Mf=min(dff)/dt*tfps_noisy.dx/tfps_noisy.dy;
                tfps_trif=apply_trifilter(tfps_noisy,Mi,Mf,1);
                a=tfps_trif.Z(:);
                a=sort(a);
                ind=round(length(a)*par.custom);
                thresh=a(ind);
                %hough e valore critico
                [peaks,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,thresh,1,1);
                [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks,Tfft,[fmin fmax],n,kn);
                pixel=Pixel(hfdf,x0,kn,gridx,gridk);
                massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
                media=mean_sub_Hough(pixel);
                dev_std=deviazione_std(pixel);
                crmat_dff(j,k,i)=(massimo-media)/dev_std;

                %filtro generico, non dipende dai parametri
                Mi=-0.01;
                Mf=fdotmax*tfps_noisy.dx/tfps_noisy.dy;
                tfps_trif=apply_trifilter(tfps_noisy,Mi,Mf,1);
                a=tfps_trif.Z(:);
                a=sort(a);
                ind=round(length(a)*par.custom);
                thresh=a(ind);
                %hough e valore critico
                [peaks,inds,times]=LAB_create_peakmap_from_wnimg(tfps_trif,thresh,1,1);
                [hfdf,gridk,gridx] = LAB_hough_and_peak_noantsour_nogd(peaks,Tfft,[fmin fmax],n,kn);
                pixel=Pixel(hfdf,x0,kn,gridx,gridk);
                massimo=massimo_intorno(hfdf,x0,kn,gridx,gridk);
                media=mean_sub_Hough(pixel);
                dev_std=deviazione_std(pixel);
                crmat_gen(j,k,i)=(massimo-media)/dev_std;

            end
        end
    end


%     %valore critico medio per strain
%     crvec=mean(crmat,[2 3]);
%     crstd=std(crmat,[],[2 3]);
%     crvec_dff=mean(crmat_dff,[2 3]);
%     crstd_dff=std(crmat_dff,[],[2 3]);
%     crvec_gen=mean(crmat_gen,[2 3]);
%     crstd_gen=std(crmat_gen,[],[2 3]);
% 
%     figure,semilogx(h0, crvec,'.','MarkerSize',10)
%     hold on
%     semilogx(h0, crvec_dff,'.','MarkerSize',10)
%     semilogx(h0, crvec_gen,'.','MarkerSize',10)
%     title('mean critical value for each h0')
%     xlabel('h0')
%     ylabel('cr')
%     ylim([min(crvec_gen(:))-1 max(crvec_gen(:))+1])
%     set(gca,'ytick',0:1:ceil(max(crvec)))
%     grid on
%     legend('normal','custom filter','generalised filter','Location','northwest')
%     hold off
% 
%     figure,errorbar(h0, crvec,crstd,'.','MarkerSize',7)
%     hold on
%     errorbar(h0, crvec_dff,crstd_dff,'.','MarkerSize',7)
%     errorbar(h0, crvec_gen,crstd_gen,'.','MarkerSize',7)
%     title('mean critical value for each h0')
%     xlabel('h0')
%     ylabel('cr')
%     ylim([min(crvec_gen(:))-1 max(crvec_gen(:))+1])
%     set(gca,'ytick',0:1:ceil(max(crvec)))
%     grid on
%     legend('normal','custom filter','generalised filter','Location','northwest')
%     hold off

    par.time=toc

    %salvo le variabili utili (potrei anche salvare crmat) e le nomina con la
    %data, per evitare sovrascritture con lo stesso nome
    nome=strcat('run_O3_',num2str(p)); 
    save(nome, 'h0','crmat','crmat_dff','crmat_gen','par');
end