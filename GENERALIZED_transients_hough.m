function [hfdf_orig, hfdf_summed, dx,gridx]=GENERALIZED_transients_hough(peaks,hm_job,braking_index,ref_perc_time,gridk,sumflag)
%HFDF_HOUGH  creates a f/df Hough map (based on 2011 version)
%
%     [hfdf,job_info,peaks,checkE]=hfdf_hough(peaks,basic_info,job_info,hm_job)
%
%
%    peaks(5,n)      peaks of the peakmap as [t,fr,amp,wnoise,wien] (corrected for the Doppler effect)
%                      time is MJD, sd in Hz/s
%    basic_info
%    job_info
%    hm_job          hough map structure 
%        .oper       'adapt' (def), 'noiseadapt' or 'noadapt'
%        .fr         sel  mode : [minf df enh maxf] min fr, enhanced frequency step, maxfr
%                    full mode or enhanced full mode: doesn't exist
%        .sd         sel mode or enhanced full mode: [minsd dsd nsd] min sd, step, number of sd
%                    full mode : doesn't exist
%        .mimaf      used for refined, otherwise nat_range is used
%        .mimaf0       "         "        
%
%    hfdf            hough map (gd2)
%    job_info        job info structure
%    checkE          service structure for test and debug

% Version 2.0 - October 2013 
% Part of Snag toolbox - Signal and Noise for Gravitational Antennas
% by Sergio Frasca - sergio.frasca@roma1.infn.it
% Department of Physics - Universit? "La Sapienza" - Rome

%tic;
% 
%peaks=job_pack_0.peaks;
%basic_info=job_pack_0.basic_info;

pow=braking_index-1;

if ~exist('ref_perc_time','var')
    ref_perc_time=0.5;
end

% 
%     if braking_index==0 || braking_index==1 %case where the sign of the constant flips
%         start_y=hm_job.kmin;
%         end_y=hm_job.kmax;
%         deltay=hm_job.kstep;
%         gridk=start_y:deltay:end_y;
%     else
%         start_y=hm_job.kmin;
%         end_y=hm_job.kmax;
%         deltay=hm_job.kstep;
%         gridk=start_y:deltay:end_y;
%     end
% end


[n1,n2]=size(peaks);
peaks(4,:)=ones(1,n2);
peaks(4,peaks(2,:)<0)=0;
%I1000=90000;%10000;  % security belt (even)
%I500=I1000/2;
Day_inSeconds=86400;
try
    epoch=quantile(peaks(1,:),ref_perc_time); %beginning, middle, or end of observing run
catch
    epoch=min(peaks(1,:)); %if the above fails, place t0 at the beginning of the observing run
end
% epoch_to_calc_I500=min(peaks(1,:));
% 
% ts=peaks(1,:);
peaks(1,:)=peaks(1,:)-epoch;


minf0=hm_job.minf;
maxf0=hm_job.maxf;
df=hm_job.df;
enh=hm_job.frenh;


if braking_index==1 % case of pulsar winds
    x=log(peaks(2,:));
    minx=min(x);
    maxx=max(x);
    dx=df*1/maxf0;
    pow=1; %not physical, just to negate contribution of pow in each expression that follows, because in exp case, there is no 'pow'
else
    x=peaks(2,:).^-pow;
    minx=1/maxf0^pow;
    maxx=1/minf0^pow; 
    
    if maxx<minx
        temp1=minx;
        temp2=maxx;
        maxx=temp1;
        minx=temp2;
    end
    
    dx=(pow)*df*1/maxf0^(braking_index); %smallest possible step
    dx=abs(dx);
end

ddx=dx/enh; %refined step in x

dx2=dx/2;

% inix=minx-dx2-ddx;
% finx=maxx+dx2+ddx;
inix=minx;
finx=maxx;
% same_ts=ts-epoch_to_calc_I500; %always using minimum time to calculate max spindown
% 
% max_time=Day_inSeconds*max(same_ts)/ddx;
% %max_time=Day_inSeconds*max_time/ddx;
% max_sd=max(yi);
% I500=ceil(1+max_sd*(max_time)*abs(pow));
% I1000=I500;%*2;

I1000=0;
I500=I1000;
% if braking_index==0 || braking_index==1 %CW case
%     I1000=4000;  % security belt (even)
%     I500=I1000/2;
% end

nbin_x=ceil((abs(finx-inix))/ddx)+I1000; %was ceil()
% deltax2=round(enh/2+0.001);
% deltax2=(dx/2)/ddx;
n_of_peaks=length(peaks);
ii=find(diff(peaks(1,:)));  % find the different times
ii=[ii n_of_peaks]; 
nTimeSteps=length(ii); % number of times of the peakmap





%yi=1e-4:1e-4:1e-2; %based on alpha=0.01-0.1
%deltay=yi(2)-yi(1);


nbin_yi=length(gridk);

ii0=1; %nbin_d,nbin_f0

nbin_yi

nbin_x

binh_df0=zeros(nbin_yi,nbin_x);  %  HM matrix container

%nbin_x=ceil((finx-inix)/ddx)+I500*2; 

%ssmax=0;
%vec=[];
for it = 1:nTimeSteps
    kf=(x(ii0:ii(it))-inix)/ddx;  % normalized xs
    
    w=peaks(4,ii0:ii(it));               % wiener weights
    t=peaks(1,ii0)*Day_inSeconds; % time conversion days to s
    tddx=t/ddx;
    x0_a=kf;
    clear kf
%     x0_a=kf-deltax2; 
%     x0_b=kf+deltax2;

    for id = 1:nbin_yi   % loop for the creation of half-differential map

        td=gridk(id)*tddx*pow;
%         a=1+round(x0_a-td+I500);  %bad security belt
%         binh_df0(id,a)=binh_df0(id,a)+w;

        inds=round(x0_a-td+I500); %good security belt
        ind_of_inds=find(inds>0);
        a=inds(ind_of_inds);
        binh_df0(id,a)=binh_df0(id,a)+w(1:length(a)); % left edge of the strips
        
        
% 
%         b_inds=round(x0_b-td+I500);
%         b_ind_of_inds=find(b_inds>=0);
%         b=1+b_inds(b_ind_of_inds);
%         binh_df0(id,b)=binh_df0(id,b)-w(1:length(b)); % right edge of the strips
%         
        
    end
    ii0=ii(it)+1;
%    ss=length(x0_a);
%     if ss>ssmax
%         ssmax=ss;
%     end 
end


hfdf_orig=gd2(binh_df0.'); 
% clear binh_df0;

hfdf_orig=edit_gd2(hfdf_orig,'dx',ddx,'ini',inix-ddx*I500,'capt','Histogram of y0-x0'); %inix-I500*ddf


if sumflag==1

    freqs_grid=minf0:df:maxf0-df; % constructs the grid in frequency
    gridx=1./freqs_grid.^(pow); % creates nonuniform grid in x with transformation
    gridx=flip(gridx); %orders the grid from smallest to largest value
    biggest_gridx=x_gd2(hfdf_orig); % overresolved grid extracted here

    % stepss=diff(gridx);
    the_rule=[];
    summed_binh_df0=[];
    for i=1:length(gridx)-1 % for each point in the grid
        first_x=gridx(i);  %x1
        second_x=gridx(i+1); %x2
        ind_first=min(find(biggest_gridx>=first_x)); % the first occurance for the closest index to the nonuni grid x1
        ind_second=max(find(biggest_gridx<=second_x)); % the last occurnce for the index closest to x2

        if i>1 %indexing stuff, so that the index used when summing corresponds to 1:3, 4:8, 9:15, 16:20 etc...
            if ind_first==the_rule(2*(i-1))
                ind_first=ind_first+1;
            end
        end

        the_rule=[the_rule ind_first ind_second]; %records the indexing rule
    %     try
        one_coln=sum(binh_df0(:,ind_first:ind_second),2); %sums the bins in x_0 according to the indices
    %     catch
    %         disp('')
    %     end
        summed_binh_df0=[summed_binh_df0 one_coln]; %creates new Hough number matrix counts
    end
    clear binh_df0;
    gridx=gridx(1:length(gridx)-1);
    hfdf_summed=gd2(summed_binh_df0.');

    clear peaks; clear x; clear summed_binh_df0

    hfdf_summed=edit_gd2(hfdf_summed,'x',gridx,'capt','Histogram of y0-x0'); %inix-I500*ddf
else
    hfdf_summed=[];
    gridx=[];
end
end